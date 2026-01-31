# Simulador de Propagación Atmosférica: Roble-Seetrue

Este repositorio contiene una herramienta escrita en Python para simular la propagación de haces de luz (Láser/LED) a través de la atmósfera y calcular la potencia recibida en un detector distante. El software considera la topografía del terreno y los efectos de la turbulencia en la geometría del haz.

### Descripción General

El programa permite visualizar y calcular cómo llega un haz de luz desde un punto A a un punto B (o con retorno mediante un retrocubo), permitiendo estudiar cómo afectan las condiciones atmosféricas y geométricas a la detección de la señal. Genera perfiles de elevación y estima la relación señal a ruido (SNR) en el sensor.

### Requisitos y Ejecución

Se requiere Python 3 y las siguientes librerías:

```bash
pip install customtkinter numpy matplotlib rasterio shapely pyproj
