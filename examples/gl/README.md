# Gibson-Low Standalone and Coupled Examples

## Standalone
Standalone runs use a spherical curvilinear cell-center grid on which to calculate the Gibson-Low solution. You can modify the dimensions/extents in the `<*dir />` sections. Time is specified in the `<CME />` section for standalone runs, `tfin` and `dt` are to be specified in seconds. 

```xml
<idir min="1.0" max="3.0" N="64"/>
<jdir min="0.01" max="0.99" N="64"/>
<kdir min="0.0" max="1.0" N="64"/>
```

- `3D64.xml` is a quick sub 1 second low resolution spheromak run to preview behvior.
- `3D256x256x64.xml` higher resolution spheromak run used to compare to SSW.
- `Tether256x256x64.xml` is an example of a tethered tear-drop geometry with collapsed points.

The output hdf5 files are in a format that can be post-processed by `genXDMF.py` to view in Paraview. 

## Coupled
Contained in `examples/helio` these runs follow similiar parameter setting as in the standalone version, but take time and grid parameters from the GAMERA set in the `<Gamera />` section. 

You can use a symmetric ```symgiblow.xml``` or use a WSA solution ```wsagiblow.xml``` as the background. 

Co-rotation of the CME is currently not-working in MPI runs as of Nov 10, 2022. This will be resolved when functionality of importing a global grid is implemented, similar to the import of the WSA solution and grid.