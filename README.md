# Extracting the full conductivity tensor in a rectangular sample

Python and Mathematica code implementing the formulas of

> Julia Gelfond and Oskar Vafek, *Extracting the full conductivity tensor in a rectangular sample*, [arXiv:2608.06643](https://arxiv.org/abs/2608.06643) (2026).

Four four-terminal resistance measurements on a rectangular sample — three with contacts at the corners and one using an edge midpoint — determine the **full 2D conductivity tensor**: the two principal conductivities σ₊ ≥ σ₋, the Hall conductivity σ_H, and the principal-axes angle α.

The original release scripts are kept unchanged at the repository root:

- `extract conductivity tensor - python.py` — standalone Python script
- `extracting conductivity tensor.nb` — equivalent Mathematica notebook

The Python code is additionally packaged as the installable module **`rectcond`** (see below), with a worked example in [`examples/extract_conductivity_tensor.ipynb`](examples/extract_conductivity_tensor.ipynb).

## Installation

```sh
pip install git+https://github.com/jag47/Extracting-the-full-conductivity-tensor-in-a-rectangular-sample.git
```

or clone and install locally (editable):

```sh
git clone https://github.com/jag47/Extracting-the-full-conductivity-tensor-in-a-rectangular-sample.git
cd Extracting-the-full-conductivity-tensor-in-a-rectangular-sample
pip install -e .
```

Dependencies (`numpy`, `scipy`, [`mpmath`](https://mpmath.org/)) are installed automatically.

## Usage

```python
from rectcond import extract_tensor

result = extract_tensor(
    d1=2.7,   # sample side along x (horizontal)
    d2=1.2,   # sample side along y (vertical)
    R1=0.1636153846153846,
    R2=0.01676923076923077,
    R3=-0.055923076923076916,
    R5=0.04815384615384615,
)

print(result.sigma_plus)   # larger principal conductivity
print(result.sigma_minus)  # smaller principal conductivity
print(result.sigma_hall)   # Hall conductivity
print(result.alpha)        # principal-axes angle, radians in [0, pi)
```

Side lengths enter only through their ratio. Conductivities are sheet conductances in units of 1/[R] (e.g. siemens per square for resistances in ohms); `alpha` is measured counterclockwise from the bottom (d1) edge.

## Measurement configurations

The sample is a rectangle with side d₁ along x and d₂ along y, corners labeled bottom-left (BL), bottom-right (BR), top-right (TR), top-left (TL). Each resistance is (Φ_A − Φ_B)/I with current I flowing from source S to drain D (Fig. 1f–h of the paper):

| Configuration | Source S | Drain D | Probe A (V+) | Probe B (V−) | Definition |
|---|---|---|---|---|---|
| 1 | BL | BR | TR | TL | R₁ = −ΔΦ₁/I (so R₁ > 0) |
| 2 | BR | TR | BL | TL | R₂ = ΔΦ₂/I |
| 3 | BL | TR | BR | TL | R₃ = ΔΦ₃/I (sign is physical; may be negative) |
| 4 (optional) | BR | TL | BL | TR | cross-check: R₄ = −2R₁ + 2R₂ − R₃ |
| 5 | BR | TR | midpoint of bottom edge | TL | R₅ = ΔΦ₅/I |

Notes:

- **R₁ is negated** so that it is conventionally positive; R₂ and R₅ come out positive as defined; **the sign of R₃ is meaningful** (it carries the Hall response) and must not be dropped.
- Configuration 4 adds no new information but provides a consistency check on the three vertex measurements (Eq. 39 of the paper).
- Consistency requires **R₂ < R₅**; if not, one of the underlying assumptions (e.g. a spatially uniform conductivity tensor) is violated.
- Configurations 1–3 determine σ_H and the geometric mean √(σ₊σ₋); the midpoint configuration 5 is what fixes the principal-axes angle α and the anisotropy ratio.

## Citing

If you use this code, please cite the paper (see [`CITATION.cff`](CITATION.cff)):

```bibtex
@misc{gelfond2026extracting,
  title  = {Extracting the full conductivity tensor in a rectangular sample},
  author = {Gelfond, Julia and Vafek, Oskar},
  year   = {2026},
  eprint = {2608.06643},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.mtrl-sci},
  url    = {https://arxiv.org/abs/2608.06643},
}
```
