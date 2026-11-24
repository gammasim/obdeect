# Python API

The package exports the optical model and the helpers below; everything else is
internal. Units are metres, nanoseconds and nanometres. For frames, see
[Optical models](optical-models.md#frames-and-units).

## Optical models

Create an optical model with `OpticalModel.from_records` from validated parameter records (pass the
air group index `propagation_group_index`; without it 1.0 is used with a DeprecationWarning), or `OpticalModel.from_simulation_models`
from a simulation-models checkout (the index comes from the site; a telescope without a site needs
it). The C++ core validates the resulting structured Python object directly. `load` reads an
optical-model JSON file and checks its `sha256`; the native loader also accepts JSON text, which
it checks the same way.

```{eval-rst}
.. autoclass:: obdeect.ray_tracing.model.OpticalModel
   :members:

.. autofunction:: obdeect.ray_tracing.model.trace_many

.. autoexception:: obdeect.ray_tracing.errors.ModelError
```

## Photon streams

```{eval-rst}
.. autofunction:: obdeect.ray_tracing.stream.trace_stream

.. autofunction:: obdeect.ray_tracing.stream.ground_to_telescope_rotation

.. autofunction:: obdeect.ray_tracing.stream.minimum_input_plane_height_m
```

## Analysis

```{eval-rst}
.. autofunction:: obdeect.ray_tracing.psf.containment

.. autofunction:: obdeect.ray_tracing.psf.summary
```
