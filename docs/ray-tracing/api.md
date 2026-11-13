# Python API

The package exports the optical model and the helpers below; everything else is
internal. Units are metres, nanoseconds and nanometres. For frames, see
[Optical models](optical-models.md#frames-and-units).

## Optical models

```{eval-rst}
.. autoclass:: obdeect.ray_tracing.model.OpticalModel
   :members:

.. autofunction:: obdeect.ray_tracing.model.trace_many

.. autoexception:: obdeect.ray_tracing.formats.ModelError
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

## Compilation

```{eval-rst}
.. autofunction:: obdeect.ray_tracing.compiler.compile_model

.. autofunction:: obdeect.ray_tracing.repository.read_telescope
```
