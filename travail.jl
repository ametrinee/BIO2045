# ---
# title: Effet du feu contrôlé sur le maintien des prairies : simulation d'un modèle de transition végétale
# repository: ametrinee/BIO2045
# auteurs:
#    - nom: Hong
#      prenom: Sumi
#      matricule: 20233311
#      github: ametrinee
#    - nom: Nguyen
#      prenom: Mathilde
#      matricule: XXXXXXXX
#      github: Mathilde389
# ---

# # Introduction

# # Présentation du modèle

# # Implémentation

# ## 1. Importation des packages nécessaires au modèle

using CairoMakie
using Distributions

import Random
Random.seed!(123456)

# ## 2. ...

(Utilise "##" pour commenter au code)

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
        timeseries[:, generation+1] .+ pop_change
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

# # Présentation des résultats

# # Discussion

# On peut aussi citer des références dans le document `references.bib`,
# @ermentrout1993cellular -- la bibliographie sera ajoutée automatiquement à la
# fin du document.

















