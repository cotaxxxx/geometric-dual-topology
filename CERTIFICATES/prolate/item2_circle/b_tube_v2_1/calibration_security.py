"""Source, workflow, namespace, and kernel security checks."""
from calibration_context import *


def _source_forbidden_code(source: str) -> list[str]:
    patterns = (
        "flo" + "at(", "Dec" + "imal(", "." + "str(",
        "arb(" + "str", "arf(" + "str", "mag(" + "str",
    )
    code = "".join(
        token.string
        for token in tokenize.generate_tokens(io.StringIO(source).readline)
        if token.type not in {tokenize.STRING, tokenize.COMMENT}
    )
    return [pattern for pattern in patterns if pattern in code]


def assert_clean_source_tree(root: Path = BTUBE_ROOT) -> None:
    offenders: dict[str, list[str]] = {}
    for path in sorted(root.rglob("*.py")):
        hits = _source_forbidden_code(path.read_text(encoding="utf-8"))
        if hits:
            offenders[path.relative_to(root).as_posix()] = hits
    if offenders:
        raise CalibrationError(f"source scan failed: {offenders}")


def assert_workflow_security(path: Path = WORKFLOW_PATH) -> None:
    text = path.read_text(encoding="utf-8")
    required = (
        "permissions:\n  contents: read", "persist-credentials: false",
        "btube-v2-1-calibration-approved-*", "github.sha", "--require-hashes",
        "--only-binary=:all:",
        "Enforce B-LOCAL binding gate before any result-bearing step",
        "binding_to_final_lambda_start",
        "B-LOCAL/B-ENTRY unpinned: workflow binding run prohibited",
    )
    if any(token not in text for token in required):
        raise CalibrationError("workflow security/authorization guard missing")
    forbidden = (
        "workflow_dispatch", "pull-requests: write", "issues: write",
        "contents: write", "persist-credentials: true",
    )
    if any(token in text for token in forbidden):
        raise CalibrationError("workflow contains forbidden write/dispatch capability")


def assert_no_stale_inputs(out_dir: Path) -> None:
    if out_dir.exists():
        raise CalibrationError("fresh-only output path already exists")
    for name in {
        "resume.json", "checkpoint.json", "calibration_records.jsonl",
        "CALIBRATION_SUMMARY.json", "DELIVERY_RECEIPT.json",
    }:
        if (HERE / name).exists():
            raise CalibrationError(f"stale calibration input present: {name}")


def assert_result_namespace(value: Any, path: str = "$") -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            if key in FORBIDDEN_RESULT_KEYS:
                raise CalibrationError(f"{path}: forbidden result key {key}")
            assert_result_namespace(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            assert_result_namespace(child, f"{path}[{index}]")
    elif isinstance(value, str) and FORBIDDEN_RESULT_PREFIX in value:
        raise CalibrationError(f"{path}: production certification string forbidden")


def _assert_repo_regular_file(path: Path, repo_root: Path = REPO_ROOT) -> Path:
    if path.is_symlink():
        raise CalibrationError("dependency path is a symlink")
    resolved_root = repo_root.resolve(strict=True)
    resolved = path.resolve(strict=True)
    try:
        resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise CalibrationError("dependency path escapes repository") from exc
    if not resolved.is_file():
        raise CalibrationError("dependency is not a regular file")
    return resolved


def load_production_kernel(repo_root: Path = REPO_ROOT):
    kernel_path = _assert_repo_regular_file(repo_root / KERNEL_RELATIVE, repo_root)
    before = sha256_hex(kernel_path.read_bytes())
    if before != KERNEL_SHA256:
        raise CalibrationError("production F/F_r kernel file-byte SHA mismatch")
    spec = importlib.util.spec_from_file_location("btube_v21_calibration_kernel", kernel_path)
    if spec is None or spec.loader is None:
        raise CalibrationError("production kernel import spec unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    after = sha256_hex(Path(module.__file__).resolve(strict=True).read_bytes())
    if after != before:
        raise CalibrationError("production kernel changed during import")
    for name in ("F_arb", "dFdr_arb"):
        function = getattr(module, name, None)
        if function is None or getattr(function, "__module__", None) != module.__name__:
            raise CalibrationError("F and F_r must be supplied by the single pinned file")
    return module, kernel_path


__all__ = [name for name in globals() if not name.startswith("__")]
