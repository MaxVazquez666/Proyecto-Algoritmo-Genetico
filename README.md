# 🧬 Algoritmo Genético Mejorado para Alineación de Secuencias de Proteínas

Este proyecto implementa y mejora un algoritmo genético (AG) orientado a la alineación múltiple de secuencias de proteínas, utilizando la matriz BLOSUM62 como función de evaluación. Se muestra cómo el algoritmo original puede perfeccionarse para obtener mejores resultados de fitness durante las generaciones evolutivas.

## 📌 Objetivo del Proyecto

El objetivo principal es **diseñar una versión mejorada del algoritmo genético original**, optimizando su capacidad de búsqueda para obtener alineaciones más eficientes y biológicamente relevantes, sin perder la integridad de las secuencias.

## 🔍 ¿Qué se mejoró?

El algoritmo original tenía limitaciones a la hora de manipular huecos (*gaps*) y residuos, lo que podía provocar:

- Pérdida o alteración del orden de los aminoácidos.
- Inconsistencia en la longitud de las secuencias.
- Mutaciones que sustituían residuos y afectaban la integridad de los datos.

Para corregir esto, mi versión mejorada introduce los siguientes cambios clave:

### ✅ Cambios implementados

1. **Cruce por máscaras de gaps**  
   Nuevo método de cruza que respeta la posición de los residuos y solo intercambia las posiciones de huecos, usando máscaras binarias:
   - `1` → Residuo
   - `0` → Gap

2. **Integridad garantizada**  
   Se validan las secuencias después de cada cruce y mutación para asegurar que los residuos originales permanecen, solo alterando los huecos.

3. **Mutación adaptativa**  
   La tasa de mutación ahora depende del número de generación, comenzando con valores altos que disminuyen conforme progresa el algoritmo.

4. **Selección por torneo**  
   Se incorporó este método para reforzar la presión selectiva hacia individuos más aptos.

5. **Gráfica comparativa Original vs Mejorado**  
   Se añade una visualización que muestra claramente la mejora en el rendimiento del algoritmo optimizado.

## ✅ Validación de Integridad

La función:
python
validar_poblacion_sin_gaps(poblacion, original_sequences)
