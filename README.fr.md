**🇫🇷 Français** · [🇬🇧 English](README.md)

# Codes de Goppa et interpolation de Lagrange

> Mémoire de Master 2 en Mathématiques Appliquées (Calcul Numérique), **Université d’Antananarivo**, 2026. Encadreur : Prof. Harinaivo ANDRIATAHINY.

Les codes de Goppa sont des codes correcteurs d’erreurs construits à partir de corps finis et de polynômes. Ils sont à la base du **cryptosystème de McEliece**, l’un des principaux candidats pour la **cryptographie post-quantique**. Ce projet les étudie, ainsi que leur lien avec l’**interpolation de Lagrange**, et propose une petite **implémentation Python faite à partir de zéro** (sans dépendance) que l’on peut exécuter et tester.

## L’idée en 30 secondes

Soit `L = {γ₁, …, γₙ}` un ensemble d’éléments distincts d’un corps fini et `g(x)` un polynôme tel que `g(γᵢ) ≠ 0`. Le **code de Goppa** Γ(L, g) est l’ensemble des mots `c` tels que

```
R_c(x) = Σ cᵢ / (x − γᵢ) ≡ 0  (mod g(x))
```

L’interpolation de Lagrange relie un mot à un polynôme : `c = (c₁,…,cₙ) → F(x) = Σ cᵢ Lᵢ(x)` avec `F(γᵢ) = cᵢ`. Pour un mot reçu `r = c + e`, le **syndrome** `S(x) = Σ rᵢ / (x − γᵢ) mod g(x)` ne dépend que de l’erreur `e`, ce qui est le point de départ des algorithmes de décodage.

Principaux résultats présentés dans le mémoire :
- dimension `k ≥ n − m·t` et distance minimale `d ≥ t + 1`, avec `t = deg g` ;
- pour les codes de Goppa **binaires irréductibles** : `d ≥ 2t + 1` (paramètres `[n, ≥ n − mt, ≥ 2t + 1]`) ;
- les codes BCH apparaissent comme cas particulier (`g(z) = z^(d−1)`).

## Ce que fait le code

`src/goppa.py` implémente, en Python pur :

| Partie | Contenu |
|--------|---------|
| Corps fini | arithmétique dans GF(2^m) (tables log/antilog) |
| Polynômes | addition, multiplication, évaluation, réduction modulo g |
| Code de Goppa | polynôme syndrome `S(x)`, matrice de contrôle `H[i][j] = γⱼⁱ / g(γⱼ)` (développée en bits), matrice génératrice (noyau de H), encodage |
| Décodage | recherche dans une table de syndromes, corrige jusqu’à `t` erreurs |
| Lagrange | interpolation sur GF(2^m) et sur les rationnels |

Exemple par défaut : **GF(8)** (construit avec `x³ + x + 1`), `L` = les 8 éléments du corps, `g(x) = x² + x + 1` (irréductible sur GF(8)). On obtient un code binaire **[8, 2, 5]** qui corrige **2 erreurs**.

> Portée : c’est une implémentation pédagogique. Le décodeur utilise une table de syndromes précalculée. Les décodeurs algébriques efficaces (**Patterson**, **Berlekamp–Massey**, Euclide étendu) sont présentés dans le mémoire comme perspectives.

## Exécuter

Nécessite Python 3.8+, aucun paquet externe.

```bash
python src/goppa.py          # démo : encodage, ajout de 2 erreurs, décodage, Lagrange
python tests/test_goppa.py   # lancer les tests (ou : python -m pytest)
```

Exemple de sortie (extrait) :

```
Binary Goppa code over GF(8): n=8, t=2, k=2, d=5
codeword c     : [1, 1, 1, 0, 1, 0, 1, 0]   S_c(x) = 0 (0 for a codeword)
received r     : [1, 0, 1, 0, 1, 0, 0, 0]   (errors at positions 1 and 6)
decoded        : [1, 1, 1, 0, 1, 0, 1, 0]   -> OK
F(L_i) == c_i for all i: True
```

## Tests

La suite de tests vérifie notamment que :
- GF(8) satisfait les axiomes de corps ;
- la matrice de contrôle binaire `H` définit **exactement le même code** que la définition `Σ cᵢ/(x − γᵢ) ≡ 0 (mod g)` (les 256 mots) ;
- `H·Gᵀ = 0`, `k ≥ n − m·t` et `d ≥ 2t + 1` ;
- **tout** motif d’erreurs de poids ≤ 2 est corrigé, pour chaque mot de code ;
- l’interpolation de Lagrange retrouve le bon polynôme, par exemple les points (0, −1), (2, 3), (4, 15) donnent `x² − 1`.

## Plan du mémoire

1. **Fondements mathématiques** : arithmétique modulaire, polynômes sur Fq, corps finis, extensions de corps.
2. **Codes correcteurs** : codes linéaires, matrices génératrice et de contrôle, code dual, codes cycliques, codes BCH.
3. **Codes de Goppa** : définition, paramètres, matrice de contrôle, codes de Goppa binaires, exemple sur F₂³, comportement asymptotique, décodage, utilisation en cryptographie.
4. **Codes de Goppa et interpolation de Lagrange** : problème d’interpolation, relation avec les codes de Goppa, exemple, avantages, limites, perspectives, algorithme et implémentation Python.

<!-- À FAIRE : après la soutenance, ajouter ici le lien du mémoire en PDF et des diapositives. -->

## Perspectives

- Décodeurs complets : **Patterson** (binaire), **Berlekamp–Massey**, **Euclide étendu**.
- Opérations plus rapides dans les corps finis.
- Applications cryptographiques (McEliece / cryptographie post-quantique).

## Structure du dépôt

```
.
├── src/goppa.py            # implémentation
├── tests/test_goppa.py     # tests
├── README.md               # version anglaise
├── README.fr.md            # version française
└── LICENSE
```

## Auteur

**Haja Tahina RABARY** — Université d’Antananarivo, Faculté des Sciences, Mention Mathématiques et Informatique.

## Licence

Publié sous [licence MIT](LICENSE).
