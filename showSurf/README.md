# showSurf

showSurf is a package that allows the visualization of vertexwise statistics on the cortical surface.

See `/ABCD_DEMOS/FEMA_showSurf_demo.m` for an end-to-end example of how to plot vertexwise data using output from the `FEMA_wrapper_DEMO.m`.

## Adding a colorbar to the display:

```matlab
    % Set colorbar
    colormap(cm);
    cb = colorbar('color', 'w');
    cb.Label.Interpreter = 'latex';
    cb.Label.String = '$z-score$';
    cb.Label.FontSize = 10;
    cb.Box = 'off';
    cb.Position = [.92 .08 .02 .8150];
    caxis(clim);
```
where:
- `cm` is a colormap (*e.g.* `cm = blueblackred();`), 
- `clim` the limits for the colorbar (*e.g.* `clim = [-fmax fmax];`).
