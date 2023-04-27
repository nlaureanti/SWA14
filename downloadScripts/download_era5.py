#!/opt/pub/spack/miniconda3/4.10.3/intel/2021.3.0/bin/python
#https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=overview
#https://cds.climate.copernicus.eu/api-how-to
import cdsapi

variables = ['mean_surface_net_long_wave_radiation_flux',
                'mean_surface_net_short_wave_radiation_flux',
                'mean_total_precipitation_rate', '10m_u_component_of_wind',
                '10m_v_component_of_wind', 'convective_snowfall_rate_water_equivalent',
                'mean_sea_level_pressure', '2m_dewpoint_temperature', '2m_temperature']

variables = [ 'mean_total_precipitation_rate' 'mean_surface_net_long_wave_radiation_flux',
                'mean_surface_net_short_wave_radiation_flux']

for year in [2001]: #range(2002,2010):

    for v in variables:
        c = cdsapi.Client()

        c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [ v
            
        ],
        'year': year,
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'day': [
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
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
    },
        f'/Volumes/A1/workdir/nicole/ERA5/Era5_{v}_{year}.nc')
