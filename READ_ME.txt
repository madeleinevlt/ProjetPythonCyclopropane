Projet de Bioinformatique Structurale
Dynamique Moléculaire du cyclopropane
  ~ DE SOUSA VIOLANTE Madeleine et KABBECH Hélène

1) Lancement de la MD dans un terminal :
- Force calculée selon une dérivée analytique :
$ python3 MD_cyclopropane_deriveesAnalytiques.py
- Force calculée selon une dérivée numérique :
$ python3 MD_cyclopropane_deriveesNumeriques.py 

(Création du fichier MD_cyclopropane.dat)

2) Visualisation des résultats sous R :
$ R
> source("plot_MD.R")

(Création des fichiers Energies_MD.png, anim_MD)
