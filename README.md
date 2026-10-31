# pmma

Macros de **ROOT** (C++) para el **tratamiento de imágenes de rayos X** de un fantoma de PMMA (acrílico), obtenidas por simulación con **GAMOS/Geant4** o en el laboratorio. Incluyen normalización, corrección de campo plano, filtro de mediana y cálculo de la relación señal/ruido (SNR).

## Objetivo

Procesar y evaluar la calidad de imágenes radiográficas: quitar el ruido y las variaciones del detector, y medir qué tan bien se distingue una estructura del fondo (SNR).

## Lenguaje y tecnologías

- **C++** como macros de [ROOT](https://root.cern/) (CERN).
- Datos en archivos `.root` con histogramas 2D (`TH2F`) llamados `RTPSPDoseHistos: Dose XY_merged`, generados por GAMOS.

## Archivos

| Archivo | Descripción |
|---|---|
| `array.C` | Lee una imagen experimental en texto (`30mAs.txt`), la pasa a un histograma 2D y la guarda como `30mAs.root`. |
| `visor.C` | Abre `30mAs_FF.root`, normaliza la imagen y la muestra en escala de grises. |
| `FFNv2_2.C` | **Corrección de campo plano** (*flat field*): combina la imagen cruda, la imagen sin objeto (OB) y la de corriente oscura (DI), y calcula la SNR. |
| `medianFilterv2_1.C` | **Filtro de mediana** de 5 × 5 para reducir el ruido. Guarda el resultado en `50.6mAs_M.root`. |
| `stripv2_1.C` | Toma una franja (*strip*) de 10 píxeles de alto, la proyecta en X y calcula la SNR. |
| `dato04.root`, `imagen_DI.root` | Datos de ejemplo. |

## Uso

Con ROOT instalado, ejecute cada macro desde esta carpeta:

```bash
root -l visor.C
```

Orden típico del procesamiento:

```
array.C ──► FFNv2_2.C (campo plano) ──► medianFilterv2_1.C (filtro) ──► stripv2_1.C (SNR)
```

## Notas

- Las macros usan nombres de archivo fijos en el código (`30mAs.root`, `OB_30mAs_60cm.root`, `50.6mAs_FF.root`, ...). La mayoría **no están en el repositorio**; solo se incluyen `dato04.root` e `imagen_DI.root`. Cambie los nombres según sus datos.
- Las macros que guardan resultados usan `TFile(..., "new")`, que **no sobrescribe** un archivo existente. Bórrelo antes de volver a ejecutar.
- Se ejecutan como macros sin nombre (todo el código entre `{ }`), por eso se corren con `root` y no se compilan con `g++`.
