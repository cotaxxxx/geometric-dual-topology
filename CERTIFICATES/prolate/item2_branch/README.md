# Item 2 — global stationary branch

長球の大域停留円枝の一意性（項目2）に向けた証明記録。

## 認証済み — 境界侵入点

`B(lambda) = F(1,lambda)` の根 `lambda_partial` について、

- `B(103/50) > 0`, `B(207/100) < 0`
- `B(206538/100000) > 0`, `B(206539/100000) < 0`
- `B'(lambda) < 0` for `lambda in [206538/100000, 206539/100000]`

したがって

    lambda_partial in (206538/100000, 206539/100000)

であり、その区間で一意。とくに `lambda_partial < 206539/100000` が
厳密に成立する。項目0の有限箱上限 `206539/100000` は境界侵入点そのもの
ではなく、それを囲う狭い有理上端である。

原証明書は `certificate_item2_boundary_entry.json`、代数導出は
`BOUNDARY_ENTRY_NFORM_NOTE.md`。

## 独立再検証 — 完了

上の5主張を、`UNVERIFIED_PROVENANCE/prolate_boundary_entry_arb.py` を
一切読まず、公開された代数形だけからゼロベースで書いた独立実装で再認証した。

- `verify_change_of_variables.py`：座標変換と測度を exact に検証
- `boundary_entry_independent.py`：4個の `B(lambda)` を Arb 包囲
- `run_enclosure.py`：4 lambda × 4 psi-band を完走し、総和だけで符号判定
- `bprime_independent.py`：Dual-over-`acb` により区間全体で `B'<0`

5主張はすべて原証明書と同符号で、包囲も相互に整合する。したがって境界侵入点の
認証結論は来歴不明ファイルに依存しない。最終機械証明書は
`independent_recheck/certificate_item2_independent.json`、実行記録は
`independent_recheck/RUN_LOG.md`。

## 未認証 — 項目2本体

- `ITEM2_FRR_SCAN_NOTE.md` — `F_rr < 0` の全面走査と証明の縮約
- `ITEM2_J_JET3_NOTE.md` — 三階jetと `J = F_rr/r` の一様負性

これらは float64 の探索結果であり、認証ではない。大域停留円枝の一意性そのもの、
すなわち `F_r` の単一符号変化は未認証である。

## ディレクトリ

- `independent_recheck/` — 現行の独立認証経路
- `exploration/` — float64 参照実装。認証経路外
- `UNVERIFIED_PROVENANCE/` — 来歴不明ファイルの隔離記録。現行認証経路外

## 残る証明義務

`J < 0`、`Q > 0` の大域化、`B < 0` の大域化、および区間演算による
`J < -7/5` の認証。境界侵入点の来歴問題は解消済みで、これらとは分離する。
