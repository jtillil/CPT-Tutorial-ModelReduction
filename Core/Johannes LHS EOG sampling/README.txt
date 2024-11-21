### README ### 

### Code(ausschnitte) für Johannes zum Artikel 
### Sample-based robust model reduction for non-linear systems biology models
### Undine Falkenhagen, Christian Himpe, Jane Knöchel, Charlotte Kloft, Wilhelm Huisinga
### https://doi.org/10.1002/pamm.202200269

Die empirical observability Gramians (EOG) werden mit den Scripts "emgr.m" und "wajima09_v6_u.m" berechnet. 

- emgr: das empirical Gramian framework von Christian Himpe
- wajima09_v6_u: enthält eine Neudefinition des Modells, damit es mit dem emgr-framework funktioniert und für Geschwindigkeit

Die Ergebnisse sind gespeichert in "results/empirical_gramians_[scenario]". 
Die Modellreduktion mit Latin hypersquare sampling (LHS) mit EOG wird gestartet mit "simulate_and_reduce_numqss.m" mit dem setting (..., 'sampling', "lhs_eog").
Wichtig zu beachten ist, dass die States in der Neudefinition nicht die gleiche Reihenfolge haben wie im ursprünglichen Modell. Deswegen wird in 'simulate_and_reduce_numqss' anhand der State- und Parameter-Namen umsortiert.
Ein Beispielaufruf, den ich für die Erstellung der Vergleichstabellen im Paper verwendet habe, steht im Script "population_error_PT_low_MAIN.m".
