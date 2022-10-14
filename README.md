==================================
 Numeriska metoder och Simulering
  Epidemimodell med simuleringar
==================================

Här har vi gjort en simulering av covid19 epidemin
med hjälp av SIR-modellen. 

I SIR-modellen tar man till hänsyn till hur många
är mottagliga, hur många som blir sjuka och hur 
många som blir resistenta/immuna.

Detta är simulerat på två sätt:

1. Genom att formulera en ODE som sedan löses med hjälp     
   av SciPy's inbyggda ODE-lösare.

2. Genom att formulera en stokastisk-matris och 
   propotionsfunktion som sedan löses med en given
   GillesPie funktion.   