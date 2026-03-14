# ---
# title: "Effet du feu contrôlé sur le maintien des prairies : simulation d'un modèle de transition végétale"
# repository: ametrinee/BIO2045
# auteurs:
#    - nom: Hong
#      prenom: Sumi
#      matricule: 20233311
#      github: ametrinee
#    - nom: Nguyen
#      prenom: Mathilde
#      matricule: 20267325
#      github: Mathilde389
# ---

# # Introduction
# Dans les prairies, les arbustes peuvent progressivement envahir les zones dominées par les graminées.
# Cette dynamique peut mener à une diminution des prairies ouvertes. Une intervention fréquemment utilisée
# en écologie est le feu contrôlé, qui affecte davantage les arbustes que les herbes (Daryanto et al., 2019). Dans ce
# projet, nous simulons cette intervention à l'aide d'un modèle de transition pour évaluer si le feu contrôlé
# permet de limiter l'expansion des arbustes et de maintenir les prairies.

# # Présentation du modèle

# # Implémentation

# ## 1. Importation des packages nécessaires au modèle

using CairoMakie
using Distributions

import Random
Random.seed!(123456)

# ## 2. ...

function check_transition_matrix!(T)
    for ligne in axes(T, 1)
        if sum(T[ligne, :]) != 1
            @warn "La somme de la ligne $(ligne) n'est pas égale à 1 et a été modifiée"
            T[ligne, :] ./= sum(T[ligne, :])
        end
    end
    return T
end

function check_function_arguments(transitions, states)
    if size(transitions, 1) != size(transitions, 2)
        throw("La matrice de transition n'est pas carrée")
    end

    if size(transitions, 1) != length(states)
        throw("Le nombre d'états ne correspond pas à la matrice de transition")
    end
    return nothing
end

function _sim_stochastic!(timeseries, transitions, generation)
    for state in axes(timeseries, 1)
        pop_change = rand(Multinomial(timeseries[state, generation], transitions[state, :]))
        timeseries[:, generation+1] .+= pop_change
    end
end

function _sim_determ!(timeseries, transitions, generation)
    pop_change = (timeseries[:, generation]' * transitions)'
    timeseries[:, generation+1] .= pop_change
end

function simulation(transitions, states; generations=500, stochastic=false)
    check_transition_matrix!(transitions)
    check_function_arguments(transitions, states)

    _data_type = stochastic ? Int64 : Float32
    timeseries = zeros(_data_type, length(states), generations + 1)
    timeseries[:, 1] = states

    _sim_function! = stochastic ? _sim_stochastic! : _sim_determ!

    for generation in Base.OneTo(generations)
       _sim_function!(timeseries, transitions, generation)
    end
    return timeseries
end

## États

## Barren, Grass, Shrubs (= Barren, Herbes, Arbustes)
s = [10, 80, 10] ## On va commencer l'état où le Grass est dominant (= 80)

states_names = ["Barren", "Grass", "Shrubs"]
states_colors = [:bisque, :lightgreen, :forestgreen]

## Scénario 1 : sans intervention
T_sans_feu = zeros(Float64, 3, 3)
## Ligne = l'état actuel / Colonne = L'état de la génration suivante
T_sans_feu[1, :] = [85, 15, 0]  ## Ici, l'état actuel = Barren, Aucun possibilité de shurbs
T_sans_feu[2, :] = [5, 80, 15]  ## Ici, l'état actuel = Grass, Grasse devient Shrubs lentement
T_sans_feu[3, :] = [2, 8, 90]   ## Ici, l'état actuel = Grass, Shrubs reste Shrubs, disparaît presque pas ou jamais

# Scénario 2 : avec intervention
T_avec_feu = zeros(Float64, 3, 3)
# Shrubs sont particulièrement vulnérables au feu
T_avec_feu[1, :] = [80, 20, 0]  ## Ici, l'état actuel = Barren, Grass pousse plus facilement après un incendie
T_avec_feu[2, :] = [5, 88, 7]   ## Ici, l'état actuel = Grass, probabilité réduite de transformation de Grass en Shrubs à cause du feu
T_avec_feu[3, :] = [15, 25, 60] ## Ici, l'état actuel = Shrubs, mais la valeur est diminue (90 -> 60) , c'est donc difficile à maintenir l'état Shrubs

# Simulations
gen = 100

det_sans_feu = simulation(T_sans_feu, s; stochastic=false, generations=gen)
det_avec_feu = simulation(T_avec_feu, s; stochastic=false, generations=gen)

## Visualisation de scénario 1 : sans intervention
f1 = Figure();
ax1 = Axis(f1[1, 1], title="Sans feu contrôlé", xlabel="Nb. générations", ylabel="Nb. parcelles")

for i in eachindex(s)
    lines!(ax1, det_sans_feu[i, :], color=states_colors[i], label=states_names[i], linewidth=4)
end

axislegend(ax1)
tightlimits!(ax1)
f1

## Visualisation de scénario 2 : avec intervention
f2 = Figure();
ax2 = Axis(f2[1, 1], title="Avec feu contrôlé", xlabel="Nb. générations", ylabel="Nb. parcelles")

for i in eachindex(s)
    lines!(ax2, det_avec_feu[i, :], color=states_colors[i], label=states_names[i], linewidth=4)
end

axislegend(ax2)
tightlimits!(ax2)
f2

# # Résultats

## **Dynamique des états sans feu contrôlé:**
# La simulation du modèle sans intervention montre une modification progressive de la distribution des parcelles
# entre les 3 états (Figure 1). Au début de la simulation, la majorité des parcelles est occupée par l’état “grass”,
# soit 80 parcelles, tandis que les états “shrubs” et “barren” comptent chacun 10 parcelles.

# Au cours des premières générations, le nombre de parcelles occupées par le “grass” diminue rapidement, passant
# d’environ 80 parcelles à environ 35 parcelles vers la 10e génération (Figure 1). Parallèlement, le nombre de
# parcelles dominées par les “shrubs” augmente fortement durant cette période initiale, passant d’environ 10 parcelles 
# à plus de 40 parcelles.

# Après environ 20 générations, les valeurs des trois états atteignent un plateau. Le nombre de parcelles dominées 
# par les “shrubs” se stabilise autour de 45-50 parcelles, tandis que les “grass” se maintiennent autour de 33 parcelles
# Les parcelles “barren” augmentent légèrement au début de la simulation, passant d’environ 10 à environ 18 parcelles, 
# puis demeurent relativement constantes pour le reste des générations simulées (Figure 1).

## **Dynamique des états avec feu contrôlé:**
# Lorsque la simulation inclut un feu contrôlé, la distribution des parcelles entre les états évolue différemment (Figure 2). 
# Au début de la simulation, la distribution initiale est identique à celle du scénario précédent, avec 80 parcelles dominées
# par les “grass” et 10 parcelles dans chacun des deux autres états.

# Durant les premières générations, le nombre de parcelles dominées par les “grass” diminue rapidement, passant d’environ 80 parcelles 
# à environ 65 parcelles vers la génération 10 (Figure 2). Après cette phase initiale, la valeur reste relativement stable autour 
# de 64-65 parcelles pour le reste de la simulation. Le nombre de parcelles dominées par les “shrubs” demeure relativement faible dans 
# ce scénario. Il augmente légèrement au début de la simulation pour atteindre environ 12 parcelles, puis diminue légèrement et se stabilise 
# autour de 11 parcelles après une vingtaine de générations (Figure 2).

# Les parcelles “barren” augmentent progressivement au début de la simulation, passant d’environ 10 parcelles à environ 24 parcelles vers 
# la génération 15. Par la suite, cette valeur demeure stable jusqu’à la fin de la simulation (Figure 2).



# # Discussion

# L’objectif de cette simulation était d’évaluer si l’introduction d’un feu contrôlé peut limiter l’expansion des arbustes (shrubs) 
# et favoriser le maintien des prairies dominées par les graminées (grass). Nous avons émis l’hypothèse que l’absence de feu 
# favoriserait la dominance des arbustes, tandis que l’introduction d’un feu contrôlé permettrait de maintenir un état dominé 
# par les graminées.

# Les résultats obtenus sans feu contrôlé montrent une diminution rapide du nombre de parcelles dominées par les graminées, 
# accompagnée d’une augmentation importante des parcelles occupées par les arbustes. Cette dynamique correspond au phénomène 
# d’encroûtement ligneux (shrub encroachment), largement documenté dans les écosystèmes de prairies, où les plantes ligneuses 
# envahissent progressivement les zones dominées par les herbacées (O’Connor et al., 2019 ; Daryanto et al., 2019). 
# Dans la simulation, cette transition conduit à une stabilisation du système vers un état dominé par les arbustes après 
# plusieurs générations. Ce résultat est cohérent avec le concept d’états stables alternatifs, selon lequel les prairies et 
# les formations arbustives peuvent représenter deux équilibres écologiques distincts. Une fois les arbustes établis, 
# ils peuvent maintenir leur dominance en exploitant plus efficacement certaines ressources, notamment l’eau et la lumière, 
# ce qui limite le développement des graminées (Daryanto et al., 2019).

# Cette transition s’accompagne également d’une diminution importante des parcelles dominées par les graminées. Dans la 
# littérature,l’expansion des arbustes est souvent associée à une réduction de la productivité herbacée et à une modification 
# de la structure des prairies (Daryanto et al., 2019). Les résultats de la simulation suivent cette tendance, puisque la 
# dominance initiale des graminées diminue fortement lorsque le feu, une perturbation écologique importante dans ces écosystèmes, 
# est absent.

# À l’inverse, l’introduction d’un feu contrôlé dans le modèle permet de maintenir une proportion importante de parcelles dominées 
# par les graminées, tandis que l’abondance des arbustes demeure relativement faible. Ce résultat concorde avec de nombreuses 
# études montrant que le feu constitue un facteur écologique majeur dans le maintien des prairies ouvertes. En perturbant davantage 
# les plantes ligneuses que les graminées, les incendies peuvent limiter l’établissement et la croissance des arbustes, empêchant 
# ainsi la conversion progressive des prairies en formations arbustives (Daryanto et al., 2019 ; Ling et al., 2023). Plusieurs 
# travaux suggèrent également que l’efficacité du feu dépend en grande partie de sa capacité à affecter les jeunes arbustes, 
# qui sont généralement plus sensibles aux perturbations que les individus adultes (Daryanto et al., 2019 ; Alfaro-Reyna et al., 2025).

# La simulation met aussi en évidence une augmentation du nombre de parcelles à sol nu (barren) lorsque le feu contrôlé est appliqué. 
# Cette tendance est également rapportée dans la littérature. Bien que le feu puisse favoriser la production de biomasse herbacée et 
# maintenir les prairies, il peut aussi entraîner une augmentation de la proportion de sol exposé, notamment immédiatement après 
# les brûlages (Daryanto et al., 2019). Cette augmentation du sol nu représente un compromis important dans la gestion des écosystèmes 
# de prairies. En effet, si le feu contribue au maintien des graminées et à la limitation des arbustes, il peut aussi réduire temporairement 
# la couverture végétale et augmenter la vulnérabilité du sol à l’érosion ou à la dégradation (Daryanto et al., 2019).

# Dans l’ensemble, les résultats de cette simulation suggèrent que le feu agit comme une perturbation écologique importante permettant de 
# maintenir la dominance des graminées et de limiter l’expansion des arbustes dans les prairies. Toutefois, ils illustrent également que 
# cette stratégie de gestion peut entraîner certains effets secondaires, notamment une augmentation du sol nu, ce qui souligne l’importance 
# d’un équilibre dans l’utilisation des feux contrôlés.

# **Limites du modèle**
# Cette simulation présente certaines limites qui doivent être considérées lors de l’interprétation des résultats. Le modèle 
# repose sur un nombre limité d’états de végétation et sur des probabilités de transition fixes, alors que les écosystèmes 
# naturels sont influencés par de nombreux facteurs supplémentaires, tels que les conditions climatiques, la disponibilité 
# en eau, le broutage par les herbivores ou la variabilité de l’intensité des incendies. De plus, le modèle ne tient 
# pas compte de la régénération à long terme des arbustes après les incendies, un phénomène observé dans plusieurs écosystèmes. 
# L’intégration de ces variables dans de futures simulations permettrait d’obtenir une représentation plus réaliste de la 
# dynamique des prairies.





# On peut aussi citer des références dans le document `references.bib`,
# @ermentrout1993cellular -- la bibliographie sera ajoutée automatiquement à la
# fin du document.

















