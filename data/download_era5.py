import cdsapi
client = cdsapi.Client()
dataset = 'reanalysis-era5-pressure-levels'
request = {
      'product_type': ['reanalysis'],
      'variable': ['geopotential',
                   'temperature', "specific_humidity",
                    "u_component_of_wind",
                    "v_component_of_wind",
                    "vertical_velocity"],
      'year': ['2024'],
      'month': ['03'],
      'day': ['01', '02'], #['01', '02', '03', '04', '05', '06', '07'],
      'time': ['00:00'], #['00:00','01:00', '02:00', '03:00','04:00','05:00', '06:00', '07:00', '08:00', '09:00','10:00', '11:00', '12:00', '13:00','14:00','15:00', '16:00', '17:00', '18:00', '19:00','20:00', '21:00', '22:00', '23:00'],
      'pressure_level': ['1000', '975', '950', '925', '900', '875', '850', '825', '800', '775', '750', '700', '650', '600', '550', '500', '450', '400', '350', '300', '250', '225', '200', '175', '150', '125', '100', '70', '50', '30', '20', '10', '7', '5', '3', '2', '1'],
      'data_format': 'netcdf', # 'grib'
}
target = 'era5_pl_sample.nc'

client.retrieve(dataset, request, target)

dataset = 'reanalysis-era5-single-levels'
request = {
      'product_type': ['reanalysis'],
      'variable': [
          "10m_u_component_of_wind",
          "10m_v_component_of_wind",
          "2m_temperature",
          "mean_sea_level_pressure",
          "toa_incident_solar_radiation",
          "total_precipitation",
                ],

      'year': ['2024'],
      'month': ['03'],
      'day': ['01', '02'],
      'time': ['00:00'],
      'data_format': 'netcdf',
}
target = 'era5_sfc_sample.nc'


client.retrieve(dataset, request, target)
