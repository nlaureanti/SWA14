#!/Volumes/A1/workdir/nicole/bin/python
###https://cds.climate.copernicus.eu/api-how-to

import cdsapi

for y in [2001,2002,2011,2012,2013,2014]:
    c = cdsapi.Client()

    c.retrieve(
    'cems-glofas-historical',
    {
        'system_version': 'version_3_1',
        'variable': 'river_discharge_in_the_last_24_hours',
        'format': 'netcdf4.zip',
        'hydrological_model': 'lisflood',
        'product_type': 'consolidated',
        'hyear': y,
        'hmonth': [
            'april', 'august', 'december',
            'february', 'january', 'july',
            'june', 'march', 'may',
            'november', 'october', 'september',
        ],
        'hday': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'area': [
            30, -170, -60,
            20,
        ],
    },
    f'Glofas_{y}.nc.zip')
