# Más allá de Máxima Verosimilitud: Introducción a la Psicometría Bayesiana
Este repositorio incluye todos los materiales impartidos en el seminario de introducción a la Psicometría Bayesiana. 
1. Scripts
   - BMPS2025_Blavaan.R: Código de R para simular datos y realizar análisis psicométricos bayesianos con `{blavaan}`. 
       1. Simulación del modelo factorial confirmatorio.
       2. Prior Predictive Chekvs
       3. Estimación del modelo
       4. Evaluación de convergencia y eficiencia
       5. Índices de ajuste total e incremental del modelo bayesiano
       6. Comparación de modelos bayesianos
       7. Resumen de resultados
       8. Posterior Predicritive Model Checks
   - BMPS2025_Stan.R: Código de R para simular datos y realizar análisis psicométricos bayesianos con `{cmdstanr}`.
       1. Simulación del modelo factorial confirmatorio
       2. Estimación del modelo
       3. Resumen de resultados
       4. Comparación de modelos
2. Slides
   - Images: imágenes utilizadas en las diapositivas de la presentación.
   - blavaan objects: objetos estiamdos con `{blavaan}` para la presentación.
   - cmdstanr objects: objetos estiamdos con `{cmdstanr}` para la presentación.
   - presentation.qmd: diapositivas de la presentación en quarto.
   - presentation.html: formato final de la presentación (no la versión pdf).
   - rrs.scss: aspectos estéticos de la presentación em html.
3. Stan models
   - CFA_marginal.stan: un ejemplo de modelo factorial confirmatorio utilizando la función de verosimilitud marginal. 
