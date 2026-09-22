"""plot a seafloor case using the installed gravmagpy plotting utilities"""

from pathlib import Path
import sys

from gravmagpy.plotting import plot_fields

if len(sys.argv) not in (3, 4):
    raise SystemExit("usage: plot_bxyz.py input.in output.txt [plot.png]")

input_path = Path(sys.argv[1]).expanduser().resolve()
output_path = Path(sys.argv[2]).expanduser().resolve()
image_path = Path(sys.argv[3]).expanduser().resolve() if len(sys.argv) == 4 else None
image_path = plot_fields(input_path, output_path, image_path)
print(f"wrote {image_path}")
