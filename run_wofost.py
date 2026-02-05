from pcse.models import Wofost72_WLP_FD
from pcse.input import (
    CSVWeatherDataProvider,
    YAMLCropDataProvider,
    YAMLAgroManagementReader
)

# Weather
weather = CSVWeatherDataProvider("weather/weather.csv")
# Crop
crop = YAMLCropDataProvider("crop")

# Soil (temporary: use None)
soil = None

# Agromanagement
agro = YAMLAgroManagementReader("agromanagement/agro.yaml")

# Model
model = Wofost72_WLP_FD(
    cropdata=crop,
    soildata=soil,
    weatherdata=weather,
    agromanagement=agro
)

model.run_till_terminate()

output = model.get_output()
final = output[-1]

print("🌾 Crop Yield (kg/ha):", final["TWSO"])
