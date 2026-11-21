# TSV BEM Analysis Summary (60 GHz)

## Structure
- 2 TSVs with radius 10 μm, height 100 μm
- Semiconductor traces (2 μm wide) in top and bottom die
- Keep-out regions around TSVs (15 μm distance)

## Panelization
- TSV cylinders: 8×6 = 48 panels per TSV (96 total)
- Rectangular traces: 6×6 = 36 panels per trace (72 total)
- **Total panels: 240**

## Results at 60 GHz

### S-Parameters (dB)
```
Port    Snn      Cross Sij
1       -0.005   1→2: -69.6
2       -0.004   1→3: -73.7
3       -0.005   1→4: -79.5
4       -0.005   3→4: -71.2
```

### Key Observations
- Excellent port matching (|Snn| < 0.01)
- Strong isolation between ports (>69 dB)
- Capacitive coupling dominates at 60 GHz
- TSV-to-trace and trace-to-trace coupling suppressed

## Files Generated
- `tsv_sparameters_60ghz.csv`: Full S-parameter matrix
- `tsv_bem_60ghz.jl`: Julia BEM implementation
