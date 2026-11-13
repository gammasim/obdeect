# Shadowing

Obscurers such as camera bodies, baffles and telescope structures block photons on their way
to the focal surface. A blocked photon records the obscurer's id in `terminal_surface`. The
obscurers are defined in [Optical elements](optical-elements.md#obscuring-elements).

## MST structure

```{figure} ../_static/figures/mst_shadowing.png
:alt: MST structure coloured by the fraction of light each element alone blocks.

MSTS-01 (FlashCam), simulation models 7.0.0: an on-axis star at infinity,
1,000,000 photons, zenith angle 20 degrees and source seed 11. [1](#figure-provenance)
```

Dotted elements are never crossed in this sample. Blue shows the fraction of unobscured light
that each element alone blocks; grey elements stop no light on their own in this sample. The
nominal layout is drawn by simtools; shadowing is calculated by obdeect ray tracing.

## Which elements matter?

`OpticalModel.obscurer_hits` traces photons without obscurers, then checks every obscurer
along paths that reach the focal surface. It reports summed input weights per element:

| Quantity | Meaning |
| --- | --- |
| `crossed` | Light whose path meets the element. |
| `first` | Light absorbed first by this element in the full optical model. |
| `sole` | Light stopped only by this element; removing it would recover this light. |

```python
from obdeect.ray_tracing import OpticalModel

model = OpticalModel.load("mst.json")
photons = model.star_photons(1_000_000, seed=11)
hits = model.obscurer_hits(photons, seed=12)
share = hits["sole"] / hits["unobscured"]
```

Command-line equivalent:

```sh
obdeect ray-tracing obscurers --model mst.json --seed 11 --photons 1000000 --field-deg 0 0
```

This diagnostic checks every obscurer on every leg and is meant for studies. An element that
blocks nothing in one sample may matter at another source position. Check the full field of
view and the internal calibration sources before removing elements. The
[telescope structure](optical-elements.md#telescope-structure) section describes the table.

## Figure provenance

[1] [MST shadowing figure](../_static/figures/telescopes-provenance.json).
