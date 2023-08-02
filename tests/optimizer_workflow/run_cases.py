from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from pyopenfoam import OpenFOAM


@dataclass
class OperationSchedule:
    flood_minutes: int
    dry_minutes: int
    end_minutes: int

    def __post_init__(self):
        self.current_time = 0.0
        self.next_time_minutes = 0.0

    @property
    def flood_dry_ratio(self):
        return Fraction(self.flood_minutes, self.dry_minutes)

    def __bool__(self):
        return self.current_time < self.end_minutes


def main():
    of = OpenFOAM("./TEMPLATE")
    schedule = OperationSchedule(flood_minutes=20, dry_minutes=60, end_minutes=500)

    print(f"{schedule = }")

    while schedule:
    # while False:

        print("\n☀️ \t Dry period", f"{schedule.current_time = }")
        print(f"{of.latest_time = }")
        
        schedule.next_time_minutes += schedule.dry_minutes
        print("Run until", f"{schedule.next_time_minutes = }")

        of.set_boundary_fixedGradient()  ## Dry
        of.set_endtime(schedule.next_time_minutes)
        of.run_case("unsatFoam")

        print("Dry period ended")
        schedule.current_time = schedule.next_time_minutes 

        if not schedule:
            break

        print("\n🌧️ \t Flood period", f"{schedule.current_time = }")
        print(f"{of.latest_time = }")

        schedule.next_time_minutes += schedule.flood_minutes
        print("Run until", f"{schedule.next_time_minutes = }")

        of.set_boundary_fixedValue()  ##Flood
        of.set_endtime(schedule.next_time_minutes)
        of.run_case("unsatFoam")

        print("Flood period ended")
        schedule.current_time = schedule.next_time_minutes

        # print("*****************")
    
    of.foam_to_vtk()

    scalar = of.read_field_all_times(field="Sw")
    
    igt = 0
    import matplotlib.pyplot as plt
    fig, (cax,ax) = plt.subplots(2,1, figsize=[5,7], gridspec_kw={"height_ratios":[0.2,5]}, sharex=False)
    img = ax.pcolormesh(scalar.t[igt:]/86400, scalar.z, scalar[:,igt:], cmap="copper")
    ax.spines.right.set_visible(False)
    ax.set_xlabel("Time $t$ [d]")
    ax.set_ylabel("Depth $z$ [m]")
    plt.colorbar(img, cax=cax, orientation="horizontal")
    # cax.set_title(fr"${label}$ {units}")
    fig.tight_layout()
    plt.savefig("field.png")

if __name__ == "__main__":
    main()
