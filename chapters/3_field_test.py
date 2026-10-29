"""Chapter 3 - Field validation on the roof.

Paper: Field validation (Table field), Methods: outdoor field test. Two twin SunPower modules in landscape, one
clean and one with the 80 mm band of one filter layer (tau_I = 0.71), recorded by an IV tracer for 16 days and
compared with the SPICE loss of the equal-height sweep (chapter 2). No LTspice.

    python chapters/3_field_test.py

Outputs: results/3_field_test/ (fig_field_test, fig_field_test_loss, fig_field_photos, summary.json = the numbers of the Table
field), results/3_field_test/data/ (pairs, daily energy, loss by irradiance and by hour)
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from soilspice import config, field, photos  # noqa: E402


def main():
    field.run(config.FIELD_CSV, config.FIELD_DIR, fig_dir=config.CH3_DIR)
    shutil.copy2(config.FIELD_DIR / "summary.json", config.CH3_DIR / "summary.json")
    if photos.PHOTO_CLEAN.exists() and photos.PHOTO_BAND.exists():
        photos.figure(config.CH3_DIR / "fig_field_photos")
    else:
        print(f"  fig_field_photos skipped: put the two photos in {photos.PHOTO_DIR}")


if __name__ == "__main__":
    main()
