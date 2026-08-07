# Item 6 Region C — finite center-cap certificate on `[1,10]`

Status: **CERTIFIED**

## Certified statement

The validated Arb/Acb computation proves

\[
A_\lambda''(v)>0
\qquad
\left(0\le v\le\frac1{20},\ 1\le\lambda\le10\right).
\]

Since \(\Psi_\lambda(0)=0\),

\[
\frac{\Psi_\lambda(w)}{w}
=
\int_0^1 A_\lambda''(tw)\,dt,
\]

and therefore

\[
\boxed{
\Psi_\lambda(w)>0
\qquad
\left(1\le\lambda\le10,\ 0<w\le\frac1{20}\right).
}
\]

## Exact block assembly

The parameter interval was divided into the nine exact adjacent blocks

\[
[1,2],[2,3],\ldots,[9,10].
\]

Every block was certified independently by exact recursive bisection of the
\((v,\lambda)\)-rectangle. The combiner checked:

- all nine block statuses are `CERTIFIED`;
- every block uses the same exact interval `v=[0,1/20]`;
- the lambda blocks are exactly adjacent from `1` to `10`;
- every block has exact rational coverage;
- the total number of terminal boxes is zero.

## Completion data

- block certificates: 9
- evaluated boxes: 6203
- certified leaves: 3106
- terminal boxes: 0
- exact rational coverage: passed on every block
- exact block adjacency: passed

The smallest certified lower endpoint occurs on

\[
v\in\left[\frac{13}{640},\frac7{320}\right],
\qquad
\lambda\in\left[\frac{153}{16},\frac{77}{8}\right],
\]

and is

\[
A_\lambda''(v)
>
0.00012893655173892051421755471961636338423925380545573.
\]

## Integrity record

- workflow run: `30175950226`
- workflow head: `c8ab0b8966f4e41cb62075e352a18960bb31af4b`
- combined artifact ID: `8624336730`
- combined artifact digest / downloaded ZIP SHA-256: `14412e60c9c007f6739859545681b531bde2e88358d9558962ec49d7aecb5242`
- combined JSON SHA-256: `f06decf0ee534a30424a3d73715af5762c7b06e854412dd2df50ffa2257b2389`
- combiner script SHA-256 used by the run: `fea3c198b81f330dd184d4390805759d890e17868d51634f9e3f689bf0d3cca1`

The full 3106-leaf evidence is preserved in the nine block artifacts. The
combined artifact records the exact block hashes, counts, adjacency checks,
and worst certified leaf.

## Scope

This closes Region C only on `1<=lambda<=10`. It does not yet certify:

- the finite center cap for `10<lambda<=100`;
- the compact interior region;
- the pole cap;
- the unbounded aspect-ratio tail;
- the full item 6 theorem.
