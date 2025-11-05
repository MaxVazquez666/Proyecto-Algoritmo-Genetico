import random
import copy
import time
import matplotlib.pyplot as plt
import blosum

blosum62 = blosum.BLOSUM(62)
NFE = 0

def get_sequences():
    seq1 = "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDRWGKLVVLGAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV"
    seq2 = "MKTLLVAAAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQKELQKQLGQKAKEL"
    seq3 = "MAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKALCVFAIN"
    return [list(seq1), list(seq2), list(seq3)]

def crear_individuo():
    return get_sequences()

def crear_poblacion_inicial(n=10):
    individuo_base = crear_individuo()
    poblacion = [[row[:] for row in individuo_base] for _ in range(n)]
    return poblacion

# --- Mutación: trabaja sobre un individuo (lista de secuencias) ---
def mutar_individuo(individuo, n_gaps, p):
    """
    Inserta hasta n_gaps por secuencia con probabilidad p.
    Devuelve un nuevo individuo (lista de listas de caracteres).
    NO se sustituyen residuos para preservar integridad (solo insertion de gaps).
    """
    nuevo_ind = []
    for sec in individuo:
        sec_copy = sec[:]  # lista de caracteres
        if random.random() < p:
            posiciones = set()
            for _ in range(n_gaps):
                pos = random.randint(0, len(sec_copy))
                while pos in posiciones:
                    pos = random.randint(0, len(sec_copy))
                posiciones.add(pos)
                sec_copy.insert(pos, '-')
        nuevo_ind.append(sec_copy)
    return nuevo_ind

def igualar_longitud_secuencias(individuo, gap='-'):
    """Recibe lista de secuencias (listas) y las rellena al largo máximo del individuo."""
    max_len = max(len(fila) for fila in individuo)
    return [fila + [gap] * (max_len - len(fila)) for fila in individuo]

def evaluar_individuo_blosum62(individuo):
    global NFE
    NFE += 1
    score = 0
    n_seqs = len(individuo)
    seq_len = len(individuo[0])
    for col in range(seq_len):
        for i in range(n_seqs):
            for j in range(i + 1, n_seqs):
                a = individuo[i][col]
                b = individuo[j][col]
                if a == '-' or b == '-':
                    score -= 4
                else:
                    score += blosum62[a][b]
    return score

def eliminar_peores(poblacion, scores, porcentaje=0.5):
    idx_ordenados = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    n_seleccionados = max(1, int(len(poblacion) * porcentaje))
    ind_seleccionados = [poblacion[i] for i in idx_ordenados[:n_seleccionados]]
    scores_seleccionados = [scores[i] for i in idx_ordenados[:n_seleccionados]]
    return ind_seleccionados, scores_seleccionados

def seleccion_torneo(poblacion, scores, tam_torneo=3):
    seleccionados = []
    for _ in range(len(poblacion)):
        torneo = random.sample(list(zip(poblacion, scores)), tam_torneo)
        seleccionado = max(torneo, key=lambda x: x[1])[0]
        seleccionados.append(copy.deepcopy(seleccionado))
    return seleccionados

def tasa_mutacion_adaptativa(gen, max_gen, min_rate=0.01, max_rate=0.15):
    return max_rate - (max_rate - min_rate) * (gen / max_gen)

# --- Cruce por máscaras (preserva residuos, solo combina posiciones de gaps) ---
def _reconstruir_desde_mask(residues, mask):
    """
    Reconstruye secuencia aplicando mask (1 -> residue, 0 -> gap).
    Ajusta la máscara si el conteo de 1s no coincide con len(residues).
    """
    n_res = len(residues)
    ones = mask.count(1)

    if ones > n_res:
        # convertir algunos 1 en 0 al azar
        diff = ones - n_res
        indices_1 = [i for i, v in enumerate(mask) if v == 1]
        for i in random.sample(indices_1, diff):
            mask[i] = 0
    elif ones < n_res:
        # convertir algunos 0 en 1 al azar
        diff = n_res - ones
        indices_0 = [i for i, v in enumerate(mask) if v == 0]
        # si no hay suficientes posiciones 0, llenamos lo posible (debería raramente ocurrir)
        if diff > len(indices_0):
            raise ValueError("No hay suficientes posiciones para ajustar la máscara.")
        for i in random.sample(indices_0, diff):
            mask[i] = 1

    out = []
    it = iter(residues)
    for m in mask:
        if m == 1:
            out.append(next(it))
        else:
            out.append('-')
    return out

def cruzar_por_mascaras_gap(seq1, seq2):
    """
    seq1, seq2: listas de caracteres (misma longitud).
    Devuelve dos nuevas listas de caracteres.
    Si las secuencias tienen longitud < 2 o son triviales, devuelve copias.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Las secuencias deben tener igual longitud para cruzar.")
    if len(seq1) < 2:
        return seq1[:], seq2[:]

    mask1 = [1 if a != '-' else 0 for a in seq1]
    mask2 = [1 if b != '-' else 0 for b in seq2]
    n = len(mask1)
    punto = random.randint(1, n - 1)

    new_mask1 = mask1[:punto] + mask2[punto:]
    new_mask2 = mask2[:punto] + mask1[punto:]

    residues1 = [a for a in seq1 if a != '-']
    residues2 = [b for b in seq2 if b != '-']

    new_s1 = _reconstruir_desde_mask(residues1, new_mask1)
    new_s2 = _reconstruir_desde_mask(residues2, new_mask2)
    return new_s1, new_s2

def cruzar_poblacion_conservada(poblacion):
    """
    Para cada par de individuos, cruza secuencia a secuencia.
    Asegura que cada par de secuencias padre esté igualado en longitud antes de cruzar.
    """
    nueva_poblacion = []
    n = len(poblacion)
    indices = list(range(n))
    random.shuffle(indices)
    for i in range(0, n - 1, 2):
        padre1 = poblacion[indices[i]]
        padre2 = poblacion[indices[i + 1]]
        hijos = []
        for s1, s2 in zip(padre1, padre2):
            # igualar longitudes entre estas dos secuencias específicas
            s1_eq, s2_eq = igualar_longitud_secuencias([s1, s2])
            # s1_eq y s2_eq son listas de caracteres
            h1, h2 = cruzar_por_mascaras_gap(s1_eq, s2_eq)
            hijos.append((h1, h2))
        hijo1 = [h[0] for h in hijos]
        hijo2 = [h[1] for h in hijos]
        nueva_poblacion.extend([hijo1, hijo2])
    if n % 2 == 1:
        # si queda impares, copiamos el último sin cruzar
        nueva_poblacion.append(copy.deepcopy(poblacion[indices[-1]]))
    return nueva_poblacion

def validar_poblacion_sin_gaps(poblacion, originales):
    for individuo in poblacion:
        for seq, seq_orig in zip(individuo, originales):
            seq_sin_gaps = [a for a in seq if a != '-']
            seq_orig_sin_gaps = [a for a in seq_orig if a != '-']
            if seq_sin_gaps != seq_orig_sin_gaps:
                return False
    return True

def obtener_best(scores, poblacion):
    if not scores:
        return None, float('-inf')
    idx_mejor = scores.index(max(scores))
    return copy.deepcopy(poblacion[idx_mejor]), scores[idx_mejor]

# =======================
# Algoritmo principal
# =======================
if __name__ == "__main__":
    generaciones = 60

    # ----- ALGORITMO ORIGINAL -----
    poblacion = crear_poblacion_inicial(10)
    original_sequences = get_sequences()
    fit_hist_ori = []
    poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]

    for gen in range(generaciones):
        scores = [evaluar_individuo_blosum62(ind) for ind in poblacion]
        if not scores:
            print("Población vacía en algoritmo original. Terminando.")
            break
        poblacion, scores = eliminar_peores(poblacion, scores)
        if not poblacion:
            print("Población eliminada completamente en original. Terminando.")
            break
        best, fitness_best = obtener_best(scores, poblacion)
        fit_hist_ori.append(fitness_best)

        # Igualar longitudes de cada individuo antes de cruzar
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
        poblacion = cruzar_poblacion_conservada(poblacion)

        # Mutar (mutar_individuo ahora espera un individuo)
        poblacion = [mutar_individuo(ind, 1, 0.12) for ind in poblacion]

        # Igualar longitudes tras la mutación
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]

    print("Validación integridad original:", validar_poblacion_sin_gaps(poblacion, original_sequences))

    # ----- ALGORITMO MEJORADO -----
    poblacion = crear_poblacion_inicial(10)
    fit_hist_mej = []
    poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]

    for gen in range(generaciones):
        scores = [evaluar_individuo_blosum62(ind) for ind in poblacion]
        if not scores:
            print("Población vacía en algoritmo mejorado. Terminando.")
            break

        poblacion = seleccion_torneo(poblacion, scores)
        best, fitness_best = obtener_best(scores, poblacion)
        fit_hist_mej.append(fitness_best)

        # Igualar longitudes antes de cruzar
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
        poblacion = cruzar_poblacion_conservada(poblacion)

        tasa = tasa_mutacion_adaptativa(gen, generaciones)
        poblacion = [mutar_individuo(ind, 1, tasa) for ind in poblacion]
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]

    print("Validación integridad mejorado:", validar_poblacion_sin_gaps(poblacion, original_sequences))

    # ----- Grafica -----
    plt.plot(fit_hist_ori, label='Original', color='blue')
    plt.plot(fit_hist_mej, label='Mejorado', color='red')
    plt.legend()
    plt.xlabel('Generación')
    plt.ylabel('Best fitness')
    plt.title('Comparación Algoritmo Genético: Original vs Mejorado')
    plt.show()
