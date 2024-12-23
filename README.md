# PSC

## Organisation des fichiers

01. cannings.ipynb : Premier notebook pour le modèle de Cannings. En fait on utilise que le Wright-Fisher, mais on a gardé le nom de Cannings pour ce notebook.

02. genealogy.html : Visualisation de la généalogie des individus dans le modèle de Cannings.

03. mutations.txt : output of reconstruct.ipynb

04. reconstruct.ipynb : Essai de reconstruction de l'arbre généalogique des individus avec l'information de mutations (ça n'a pas marché).

05. nca.ipynb : Notebook pour l'étude de la dynamique des mutations, notamment la statistique d'alleles non fixés (erronément appelée NCA - non coalesced alleles). 

06. nca.png : Figure de la statistique NCA.

07. wright_fisher.ipynb : Notebook pour l'étude du modèle de Wright-Fisher, qui plot encore la statistique NCA.

08. coalesce.py : fichier trouvé sur internet pour faire la coalescence des arbres généalogiques.

09. coalesce.ipynb : Notebook pour essayer de faire la coalescence manuellement, avec mutations.

10. coalesce2.ipynb : Coalescence + non fixated mutations, encore avec mutations.

11. msprime.ipynb : Notebook pour essayer de faire la coalescence avec msprime. Ça marche mieux que les autres méthodes.

12-15. wf.ipynb à wf4.ipynb : Notebooks brouillons pour le modèle de Wright-Fisher. Sur wf 4 on commence à faire des maths.

16. msprime2.ipynb : Notebook pour essayer de faire la coalescence avec msprime, mais sans mutations. Nous regardons seulement le nombre et taille des branches.
 
17. output.png : A beautiful graph of the coalescent expected value as a function of time.

18. coalescent.ipynb : On a réussi a calculer la loi hypoexponentielle que donne l'espérance du modèle de coalescence. On plot quelques graphiques sympas.

19. wf5.ipynb : document vide ??

20. lwf: on essaye de faire un modèle de Wright-Fisher à main mais Lazy, qui ne calcule pas tout à chaque itération, mais seulement quand c'est nécessaire.

## Offre PSC

Offre PSC - Mathématiques Appliquées (MAP)
Service : CMAP
École Polytechnique
Point de contact RH : Grégoire Allaire
Téléphone : 01 69 33 46 11
Site web : https://cmap.ip-paris.fr/
E-mail : gregoire.allaire@polytechnique.fr

Titre du projet : Modèles de méta-populations et résistance aux antibiotiques

Description du sujet / objectifs :
Encadrants : Vincent Bansaye et Gaël Raoul (CMAP, École Polytechnique).

La résistance aux antibiotiques est un problème de santé majeur, dans un contexte où très peu de nouvelles classes d’antibiotiques peuvent être espérées. Ces bactéries sont typiquement des bactéries commensales : ce sont des bactéries qui sont présentes à nos côtés au quotidien, en particulier dans nos intestins. Ces bactéries traversent parfois les parois de l’intestin, menant à des infections. Ces infections sont dangereuses et, lorsque les bactéries impliquées sont résistantes aux antibiotiques, elles peuvent être fatales.

Ces bactéries peuvent également causer des infections nosocomiales : profitant de la faiblesse du système immunitaire des patients ou de procédures médicales invasives, elles infectent des personnes hospitalisées. Ces infections sont particulièrement dangereuses lorsqu’elles s’établissent dans l’environnement hospitalier (équipements médicaux, par exemple) et que les bactéries impliquées possèdent de fortes capacités de résistance aux antibiotiques.

Le contrôle de la résistance des bactéries apparaît comme une problématique sociétale importante. De nombreux travaux existent pour comprendre comment une résistance aux antibiotiques peut émerger ou comment optimiser l’utilisation d’antibiotiques pour traiter des patients. Ces travaux ont des applications concrètes en termes de recommandations pour l’utilisation des médicaments ou de protocoles de soins. Malgré ces efforts, les bactéries résistantes restent très présentes et de nouveaux variants multi-résistants continuent d’émerger.

Plus récemment, une approche plus globale du problème est apparue : le traitement d’un patient n’est pas la seule échelle à considérer pour comprendre ce problème. Les échanges de bactéries entre différents lieux ainsi qu’entre humains et animaux doivent être pris en compte pour comprendre la prévalence des bactéries résistantes. Les bactéries (ou des séquences d’ADN de celles-ci) sont couramment échangées entre hôpitaux, communautés (populations humaines hors hôpitaux) et populations animales.

Ce point de vue est désigné par l’appellation One Health (Hernando-Amado S, Coque TM, Baquero F, Martinez JL. Defining and combating antibiotic resistance from One Health and Global Health perspectives. Nature Microbiology. 2019 Sep;4(9):1432-42.) et fait l’objet de projets de collecte de données importants. Les premiers résultats de ces programmes montrent que des bactéries résistantes sont présentes dans de nombreux environnements.

Quelles actions pourraient permettre de faire reculer le taux de bactéries résistantes ?
Jusqu’à présent, très peu de modèles mathématiques ont été proposés pour accompagner les biologistes et médecins travaillant sur les problématiques One Health (Niewiadomska AM, Jayabalasingham B, Seidman JC, Willem L, Grenfell B, Spiro D, Viboud C. Population-level mathematical modeling of antimicrobial resistance: a systematic review. BMC Medicine. 2019 Dec;17:1-20).

Nous pensons que des modèles simples pourraient apporter un point de vue intéressant sur ce problème.

Complexité du problème :
Différentes échelles : de la dynamique des bactéries dans un intestin aux échanges internationaux de bactéries.
Différents compartiments : hôpitaux, communautés, élevages animaux.
Différents traitements : l’utilisation d’antibiotiques varie d’un compartiment à l’autre, tout comme les conditions d’hygiène.
Échanges : les bactéries sont transportées par les individus, les aliments ou les eaux usées.
Nous proposons de considérer une partie de ces propriétés et de dériver des modèles mathématiques adaptés. Typiquement, on pourra imaginer un modèle structurant la population en trois compartiments (hôpitaux, communautés, élevages animaux) avec un traitement antibiotique spécifique à chaque compartiment, tandis que les autres propriétés décrites ci-dessus seraient négligées ou réduites à une expression très simple dans le modèle.

On s’attend alors à obtenir un modèle de méta-population. La prévalence de résistance peut alors être décrite par un système d’équations différentielles ou d’équations aux dérivées partielles (Ducrot A, Griette Q, Liu Z, Magal P. Differential Equations and Population Dynamics I. Springer, Cham: 2022).

Un travail de modélisation permettra d’écrire ce modèle, qui pourra ensuite être simulé numériquement et analysé théoriquement : existence et stabilité des points d’équilibre, influence des différents paramètres, etc. Une fois ce modèle déterministe établi, une approche stochastique pourrait être imaginée afin de suivre la dynamique des pathogènes résistants d’un compartiment à l’autre (Montagnon P. A stochastic SIR model on a graph with epidemiological and population dynamics occurring over the same time scale. Journal of Mathematical Biology, 2019 Jul;79:31-62).

Là aussi, des approches de simulation peuvent être imaginées, et des méthodes d’analyse, telles que la théorie des grandes déviations, pourraient permettre de construire une compréhension des histoires de vie des pathogènes.

Nom et Prénom : Gaël Raoul
Téléphone : 01 69 33 45 66
E-mail : gael.raoul@polytechnique.edu

Ressources nécessaires pour la réalisation du projet :
Ressources / budget mis à disposition par l’entreprise : À définir.
Estimation des risques : Aucun risque identifié.
Confidentialité : Non.
Sécurité défense : Non.
Explications : Travail théorique.