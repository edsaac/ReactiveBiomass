from dataclasses import dataclass
from fractions import Fraction

import pyopenfoam as of

of.PATH_TO_CASE = "TEMPLATE"


@dataclass
class OperationSchedule:
    flood_minutes: int
    dry_minutes: int
    end_minutes: int

    def __post_init__(self):
        self.current_time = 0.0
        self.next_time_minutes = self.dry_minutes

    @property
    def flood_dry_ratio(self):
        return Fraction(self.flood_minutes, self.dry_minutes)

    def __bool__(self):
        return self.current_time <= self.end_minutes

    def __str__(self) -> str:
        return f"""
        {self.flood_dry_ratio = }
        {self.flood_minutes = }
        {self.dry_minutes = }
        """


def main():
    schedule = OperationSchedule(flood_minutes=40, dry_minutes=100, end_minutes=1500)

    print(schedule)

    while schedule:
        print("Dry", schedule.current_time)
        schedule.next_time_minutes += schedule.dry_minutes
        print("Run until", schedule.next_time_minutes)

        of.set_boundary_fixedGradient()  ## Dry
        of.set_endtime(schedule.next_time_minutes)
        of.run_case("unsatFoam")

        schedule.current_time = schedule.next_time_minutes

        if not schedule:
            break

        print("Flood", schedule.current_time)
        schedule.next_time_minutes += schedule.flood_minutes
        print("Run until", schedule.next_time_minutes)

        of.set_boundary_fixedValue()  ##Flood
        of.set_endtime(schedule.next_time_minutes)
        of.run_case("unsatFoam")

        schedule.current_time = schedule.next_time_minutes

        print("*****************")


if __name__ == "__main__":
    main()
