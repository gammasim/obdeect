# obdeect

Installs all obdeect components at the same version:

| Component | Package | Command |
|---|---|---|
| photon stream | `obdeect-stream` | `obdeect stream ...` |
| atmospheric propagation | `obdeect-propagation` | `obdeect propagation ...` |
| light sources | `obdeect-sources` | `obdeect sources ...` |
| ray tracing | `obdeect-ray-tracing` | `obdeect ray-tracing ...` |

```bash
pip install obdeect             # everything
pip install obdeect-sources     # one component (with obdeect-stream)
```
