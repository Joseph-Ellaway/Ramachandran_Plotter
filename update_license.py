import re
import datetime

def update_license_date(license_file="./LICENSE"):
    """
    Update the year in the license file to the current year.
    """

    current_year = datetime.datetime.now().year

    # Get old year
    with open(license_file, 'r') as file:
        data = file.readlines()

        for line in data:
            # Use regex to find 202*
            match = re.search(r'20\d{2}', line)
            previous_year = match.group() if match else None
            print(previous_year)

    # No change
    if previous_year == str(current_year):
        return

    # Update the year
    with open(license_file, 'w') as file:
        for line in data:
            if previous_year in line:
                file.write(line.replace(previous_year, str(current_year)))
            else:
                file.write(line)

if __name__ == '__main__':
    update_license_date()
