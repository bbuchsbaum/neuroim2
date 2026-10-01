# Decoded accessor parity and storage-aware iteration

Epic: `bd-01M3BXX4RSBE1R4MATR8BPRGG9`. Revised 2026-09-25 after user review.

[Contract v2](block-streaming-design.md) supersedes the earlier read_block/execution-engine design. Mote owns status and acceptance; this file is a plan index.

## Now

Fix mapped decoding through existing linear_access, document the generic contract, and enforce the focused backend × datatype × public-accessor matrix. Existing scaled-integer I/O fixtures and lazy group/time operations remain the foundation.

## Then

Measure storage direction, add one deflist-style storage-aware vec_blocks iterator, qualify its file/sequence paths, and document caller-owned computation. No map/reduce engine, searchlight executor, output writer or parallel infrastructure.

## Revised tasks

| Task | Mote ID | Prerequisites |
| --- | --- | --- |
| BS-01: Correct the design: linear_access parity and storage-aware iteration | `bd-01M3BXXD6VXRKRVZ5GY5BA7P6A` | None |
| BS-08: Fix mapped physical-value decoding through linear_access | `bd-01M3BXXDRZGCY1XPE0JTDDWVX1` | BS-01 |
| BS-11: Audit existing accessor contracts across 4D and 5D representations | `bd-01M3BXXDZSTHBD9H66859C5ZKB` | BS-01 |
| BS-02: Enforce backend x datatype x public-accessor parity | `bd-01M3BXXD954YQ3AHXTEQ5F3GDM` | BS-08, BS-11 |
| BS-03: Measure backend read direction and source preparation costs | `bd-01M3BXXDC9QNETH345BZCRYVSM` | BS-01 |
| BS-13: Implement a lazy storage-aware vec_blocks iterator | `bd-01M3BXXE5ZED5C9QT8ZVEK7JV9` | BS-03, BS-02 |
| BS-09: Add volume-oriented mapped and file-backed iterator paths | `bd-01M3BXXDV3S0CSY6P9VYGEFRAM` | BS-13 |
| BS-10: Compose lazy iteration across sequences without eager concatenation | `bd-01M3BXXDXEP3QBQV8GAVQ3FQ9Q` | BS-09 |
| BS-25: Qualify iterator correctness, laziness, I/O and working memory | `bd-01M3BXXFFP8GCSRMATHE0NDB7V` | BS-10, BS-02 |
| BS-26: Document caller-owned iteration and the corrected storage boundary | `bd-01M3BXXFMRYFX232QD6KGX3VVE` | BS-25 |
| BS-27: Run package and installed-artifact gates for read parity and iteration | `bd-01M3BXXFQ84EA844G470P7D7MZ` | BS-26 |

## Withdrawn version-1 work

These tickets remain in Mote as closed cancellations with original history. They are not implemented, and independent constructor/write defects are not declared fixed.

BS-04, BS-05, BS-06, BS-07, BS-12, BS-14, BS-15, BS-16, BS-17, BS-18, BS-19, BS-20, BS-21, BS-22, BS-23, BS-24.

The replacement graph uses blocking edges only for real revised prerequisites.
Historical containment relations remain. BS-01, BS-08, BS-11, BS-02, BS-03 and
BS-13 have been implemented and locally checked. The shipped iterator supports
concrete dense, sparse and BigNeuroVec sources. Mapped/file-backed adapters
(BS-09), sequence composition (BS-10), and their subsequent qualification and
documentation gates remain open. Mote holds the detailed evidence and status.

```sh
mote show bd-01M3BXX4RSBE1R4MATR8BPRGG9
mote ls --tag scope-v2
mote ready
```

The 0.20.0 audit publishes the completed accessor and iterator subset under the
user's explicit commit/push instruction. It does not close the remaining epic
gates. Preserve concurrent work and reserve exact files before implementation.
