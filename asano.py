import math
import heapq
from bisect import bisect_left
from sortedcontainers import SortedList, SortedKeyList
from collections import deque
from exceptions import InvalidPolygonError
import ast

import random
import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

EDGE_END   = 0
EDGE_START = 1
REFLEX     = 2

def validate(v: list) -> None:
    """
    Перевірка введених даних на правильність.
    """
    n = len(v)
    for i in range(n):
        vertex = v[i]
        i_next = (i + 1) % n
        vertex_next = v[i_next]
        if vertex[0] != vertex_next[0] and vertex[1] != vertex_next[1]:
            raise InvalidPolygonError("Багатокутник містить діагональні ребра")


def fix_orientation(v: list, is_hole=False) -> None:
    """
    Фіксуємо канонічну орієнтацію багатокутника
    Якщо це зовнішня оболонка, то орієнтація буде
        проти годинникової стрілки.
    Якшо це отвір, то за годинниковою стрілкою
    """
    A = 0 # A > 0 значить, що орієнтація проти годинникової стрілки
    n = len(v)
    for i in range(n):
        vertex = v[i]
        i_next = (i + 1) % n
        vertex_next = v[i_next]
        A += vertex[0]*vertex_next[1]-vertex[1]*vertex_next[0]
    if is_hole and A > 0:
        v.reverse()
    if not is_hole and A < 0:
        v.reverse()


def determine_convexity(v: list, is_hole: bool = False) -> None:
    """
    After this call, each vertex is [x, y, is_reflex_bool],
    where is_reflex_bool == True means a reflex vertex with respect
    to the interior of the polygon (including holes).
    """
    n = len(v)
    for i in range(n):
        i_prev = (i - 1) % n
        i_next = (i + 1) % n

        vx_prev = v[i_prev]
        vx_curr = v[i]
        vx_next = v[i_next]

        # Take only first two coordinates, ignore any extras
        x_prev, y_prev = vx_prev[0], vx_prev[1]
        x_cur,  y_cur  = vx_curr[0], vx_curr[1]
        x_next, y_next = vx_next[0], vx_next[1]

        prev_curr = (x_cur - x_prev, y_cur - y_prev)
        curr_next = (x_next - x_cur, y_next - y_cur)

        turn = prev_curr[0] * curr_next[1] - prev_curr[1] * curr_next[0]

        if is_hole:
            # Hole is oriented CW by fix_orientation.
            # Convex for the hole (turn < 0) => reflex for the polygon.
            is_reflex = (turn < 0)
        else:
            # Outer is oriented CCW by fix_orientation.
            # Reflex corner for the polygon when turn < 0.
            is_reflex = (turn < 0)

        # Overwrite the vertex in a canonical form
        v[i] = [x_cur, y_cur, is_reflex]
   

    
def classify_edges(outer: list, holes: list) -> (list, list, list):
    """
    Класифікація ребер багатокутника: горизонтальні чи вертикальні.
    Збираємо всі вігнуті вершини.

    Формат вершин у outer та holes:
        vertex = [x, y, is_reflex_bool]

    Повертає:
        vertical_edges: список словників
            {
                "id": int,                 # глобальний ID вертикального ребра
                "x": float,
                "y_low": float,
                "y_high": float,
                "loop_index": int,         # 0 для outer, 1.. для дірок
                "is_outer": bool,          # True лише для зовнішньої межі
                "i": int,                  # індекс початкової вершини в цьому контурі
                "upward": bool,            # чи йде ребро вгору (від (i) до (i_next))
            }

        horizontal_edges: список словників
            {
                "id": int,
                "y": float,
                "x_low": float,
                "x_high": float,
                "loop_index": int,
                "is_outer": bool,
                "i": int,
                "left_to_right": bool,     # напрямок по x
            }

        reflex_vertices: список словників
            {
                "id": int,
                "x": float,
                "y": float,
                "loop_index": int,
                "is_outer": bool,
                "i": int,                  # індекс вершини в контурі
            }
    """

    vertical_edges = []
    horizontal_edges = []
    reflex_vertices = []

    def process_loop(vertices: list, loop_index: int, is_outer: bool):
        n = len(vertices)
        for i in range(n):
            vx = vertices[i]
            vx_next = vertices[(i + 1) % n]

            x, y, is_reflex = vx
            x_next, y_next, _ = vx_next

            # Записуємо вігнуті вершини
            if is_reflex:
                vtx = {
                    "x": x,
                    "y": y,
                    "loop_index": loop_index,
                    "is_outer": is_outer,
                    "i": i,
                }
                vtx["id"] = len(reflex_vertices)
                reflex_vertices.append(vtx)

            # Класифікація ребра
            if x == x_next:  # вертикальне ребро
                y_low = min(y, y_next)
                y_high = max(y, y_next)
                upward = (y_next > y)
                edge = {
                    "x": x,
                    "y_low": y_low,
                    "y_high": y_high,
                    "loop_index": loop_index,
                    "is_outer": is_outer,
                    "i": i,
                    "upward": upward,
                }
                edge["id"] = len(vertical_edges)
                vertical_edges.append(edge)
            else:            # горизонтальне ребро
                x_low = min(x, x_next)
                x_high = max(x, x_next)
                left_to_right = (x_next > x)
                edge = {
                    "y": y,  # y однакове для обох кінців
                    "x_low": x_low,
                    "x_high": x_high,
                    "loop_index": loop_index,
                    "is_outer": is_outer,
                    "i": i,
                    "left_to_right": left_to_right,
                }
                edge["id"] = len(horizontal_edges)
                horizontal_edges.append(edge)

    # Зовнішній контур має loop_index = 0
    process_loop(outer, loop_index=0, is_outer=True)

    # Отвори мають loop_index = 1, 2, ...
    for h_idx, hole in enumerate(holes, start=1):
        process_loop(hole, loop_index=h_idx, is_outer=False)

    return vertical_edges, horizontal_edges, reflex_vertices



# ----------------------------------------------------------------------
# 1. Event builders (separate for x / y)
# ----------------------------------------------------------------------

def build_events_y(vertical_edges: list[dict], reflex_vertices: list[dict]):
    """
    Події для горизонтального сканування (по осі y).

    Типи подій:
        EDGE_END   (0) - виходимо з вертикального ребра
        EDGE_START (1) - входимо в вертикальне ребро
        REFLEX     (2) - вігнута вершина
    """
    events = []

    for e in vertical_edges:
        # ребро активно на [y_low, y_high)
        events.append((e["y_high"], EDGE_END,   e["id"]))
        events.append((e["y_low"],  EDGE_START, e["id"]))

    for v in reflex_vertices:
        events.append((v["y"], REFLEX, v["id"]))

    events.sort(key=lambda ev: (ev[0], ev[1]))
    return events


def build_events_x(horizontal_edges: list[dict], reflex_vertices: list[dict]):
    """
    Події для вертикального сканування (по осі x).

    Типи подій:
        EDGE_END   (0) - виходимо з горизонтального ребра
        EDGE_START (1) - входимо в горизонтальне ребро
        REFLEX     (2) - вігнута вершина
    """
    events = []

    for e in horizontal_edges:
        # ребро активно на [x_low, x_high)
        events.append((e["x_high"], EDGE_END,   e["id"]))
        events.append((e["x_low"],  EDGE_START, e["id"]))

    for v in reflex_vertices:
        events.append((v["x"], REFLEX, v["id"]))

    events.sort(key=lambda ev: (ev[0], ev[1]))
    return events



# ----------------------------------------------------------------------
# 2. Horizontal / vertical sweeps: visibility "chords"
# ----------------------------------------------------------------------

def horizontal_sweep(vertical_edges: list[dict],
                     reflex_vertices: list[dict]) -> list[dict]:
    """
    Будує горизонтальні відрізки видимості (\"хорди\") за допомогою
    горизонтального сканування по y.

    На кожній висоті y підтримуємо множину активних вертикальних ребер,
    відсортовану за x. В середині інтервалів між ними записуємо, де ми
    \"всередині\" багатокутника. Для кожної вігнутої вершини v беремо
    найближчі вертикальні межі зліва/справа, які утворюють внутрішній
    інтервал, і створюємо по одному або два горизонтальні відрізки
    видимості.

    Повертає список сегментів:
        seg = {
            "id": int,
            "y": float,
            "x_from": float,
            "x_to": float,
            "vertex_id": int,
            "direction": "left" | "right",
        }
    """
    events = build_events_y(vertical_edges, reflex_vertices)

    # Активні вертикальні ребра, відсортовані за x
    active_edges = SortedKeyList(key=lambda e: e["x"])
    
   
def horizontal_sweep(vertical_edges: list, reflex_vertices: list) -> list[dict]:
    events = build_events_y(vertical_edges, reflex_vertices)

    active_edges = SortedKeyList(key=lambda e: e["x"])
    edge_by_id = {e["id"]: e for e in vertical_edges}
    vtx_by_id = {v["id"]: v for v in reflex_vertices}

    visibility_segments: list[dict] = []

    for y, kind, obj_id in events:
        if kind == EDGE_START:
            # ребро стає активним
            active_edges.add(edge_by_id[obj_id])

        elif kind == EDGE_END:
            # ребро втрачає активність
            active_edges.remove(edge_by_id[obj_id])

        else:  # REFLEX
            if not active_edges:
                continue

            v = vtx_by_id[obj_id]
            x_v = v["x"]

            xs = [e["x"] for e in active_edges]
            m = len(xs)
            if m < 2:
                continue

            inside = [(i % 2 == 0) for i in range(m - 1)]
            j = bisect_left(xs, x_v)

            if j < m and xs[j] == x_v:
                if j > 0 and inside[j - 1]:
                    visibility_segments.append({
                        "id": len(visibility_segments),
                        "y": y,
                        "x_from": x_v,
                        "x_to": xs[j - 1],
                        "vertex_id": v["id"],
                        "direction": "left",
                    })
                if j < m - 1 and inside[j]:
                    visibility_segments.append({
                        "id": len(visibility_segments),
                        "y": y,
                        "x_from": x_v,
                        "x_to": xs[j + 1],
                        "vertex_id": v["id"],
                        "direction": "right",
                    })
            else:
                if 0 < j < m:
                    idx = j - 1
                    if inside[idx]:
                        visibility_segments.append({
                            "id": len(visibility_segments),
                            "y": y,
                            "x_from": x_v,
                            "x_to": xs[idx],
                            "vertex_id": v["id"],
                            "direction": "left",
                        })
                        visibility_segments.append({
                            "id": len(visibility_segments),
                            "y": y,
                            "x_from": x_v,
                            "x_to": xs[idx + 1],
                            "vertex_id": v["id"],
                            "direction": "right",
                        })

    return visibility_segments


def vertical_sweep(horizontal_edges: list, reflex_vertices: list) -> list[dict]:
    events = build_events_x(horizontal_edges, reflex_vertices)

    active_y = SortedList()
    edge_by_id = {e["id"]: e for e in horizontal_edges}
    vtx_by_id = {v["id"]: v for v in reflex_vertices}

    visibility_segments: list[dict] = []

    for x, kind, obj_id in events:
        if kind == EDGE_START:
            e = edge_by_id[obj_id]
            active_y.add(e["y"])

        elif kind == EDGE_END:
            e = edge_by_id[obj_id]
            active_y.remove(e["y"])

        else:  # REFLEX
            if not active_y:
                continue

            v = vtx_by_id[obj_id]
            y_v = v["y"]

            ys = list(active_y)
            m = len(ys)
            if m < 2:
                continue

            inside = [(i % 2 == 0) for i in range(m - 1)]
            j = active_y.bisect_left(y_v)

            if j < m and ys[j] == y_v:
                if j > 0 and inside[j - 1]:
                    visibility_segments.append({
                        "id": len(visibility_segments),
                        "x": x,
                        "y_from": y_v,
                        "y_to": ys[j - 1],
                        "vertex_id": v["id"],
                        "direction": "down",
                    })
                if j < m - 1 and inside[j]:
                    visibility_segments.append({
                        "id": len(visibility_segments),
                        "x": x,
                        "y_from": y_v,
                        "y_to": ys[j + 1],
                        "vertex_id": v["id"],
                        "direction": "up",
                    })
            else:
                if 0 < j < m:
                    idx = j - 1
                    if inside[idx]:
                        visibility_segments.append({
                            "id": len(visibility_segments),
                            "x": x,
                            "y_from": y_v,
                            "y_to": ys[idx],
                            "vertex_id": v["id"],
                            "direction": "down",
                        })
                        visibility_segments.append({
                            "id": len(visibility_segments),
                            "x": x,
                            "y_from": y_v,
                            "y_to": ys[idx + 1],
                            "vertex_id": v["id"],
                            "direction": "up",
                        })

    return visibility_segments
# ----------------------------------------------------------------------
# 3. Побудова перетинного дводольного графа (H vs V)
# ----------------------------------------------------------------------

def build_bipartite_graph(horizontal_segs: list[dict],
                          vertical_segs: list[dict]):
    """
    Створює дводольний граф перетинів між горизонтальними та вертикальними
    відрізками видимості.

    Вхід:
        horizontal_segs: список словників
            { "y", "x_from", "x_to", ... }

        vertical_segs: список словників
            { "x", "y_from", "y_to", ... }

    Вихід:
        adj_H: список суміжності для H (горизонтальні сегменти)
               adj_H[h_id] = список v_id, з якими перетинається h_id
        adj_V: список суміжності для V (вертикальні сегменти)
    """
    H = len(horizontal_segs)
    V = len(vertical_segs)

    # Нормалізація: гарантуємо, що x_from <= x_to, y_from <= y_to
    for h in horizontal_segs:
        x1, x2 = h["x_from"], h["x_to"]
        if x2 < x1:
            x1, x2 = x2, x1
        h["_x_low"] = x1
        h["_x_high"] = x2

    for v in vertical_segs:
        y1, y2 = v["y_from"], v["y_to"]
        if y2 < y1:
            y1, y2 = y2, y1
        v["_y_low"] = y1
        v["_y_high"] = y2

    # Події по x:
    #   0 = H_START, 1 = V_QUERY, 2 = H_END
    events = []

    for h_id, h in enumerate(horizontal_segs):
        events.append((h["_x_low"], 0, h_id))
        events.append((h["_x_high"], 2, h_id))

    # ВАЖЛИВО: вертикальні сегменти мають бути відсортовані за x,
    # щоб граф став Y-конвексний на стороні V (для Imai–Asano).
    vertical_segs_sorted = sorted(
        enumerate(vertical_segs),
        key=lambda kv: kv[1]["x"]
    )
    # Перенумеровуємо v_id відповідно до цього порядку
    id_map = {}
    for new_id, (old_id, seg) in enumerate(vertical_segs_sorted):
        id_map[old_id] = new_id

    for old_id, seg in vertical_segs_sorted:
        x = seg["x"]
        v_id = id_map[old_id]
        events.append((x, 1, v_id))

    events.sort(key=lambda ev: (ev[0], ev[1]))

    # Активні горизонтальні сегменти, відсортовані за y
    # Елемент = (y, h_id)
    active = SortedList()

    adj_H = [[] for _ in range(H)]
    adj_V = [[] for _ in range(V)]

    for x, kind, obj_id in events:
        if kind == 0:
            # H_START
            h_id = obj_id
            y = horizontal_segs[h_id]["y"]
            active.add((y, h_id))

        elif kind == 2:
            # H_END
            h_id = obj_id
            y = horizontal_segs[h_id]["y"]
            active.remove((y, h_id))

        else:
            # V_QUERY
            v_id = obj_id
            if not active:
                continue

            # знайти оригінальний сегмент за позицією в відсортованому масиві
            old_id, seg = vertical_segs_sorted[v_id]
            v = seg
            y_low = v["_y_low"]
            y_high = v["_y_high"]

            # Знаходимо всі активні горизонтальні з y в [y_low, y_high]
            lo = active.bisect_left((y_low, -math.inf))
            hi = active.bisect_right((y_high, math.inf))

            for idx in range(lo, hi):
                y, h_id = active[idx]
                adj_H[h_id].append(v_id)
                adj_V[v_id].append(h_id)

    return adj_H, adj_V


# ----------------------------------------------------------------------
# 4. Стискання до конвексного представлення (інтервали) + matching
# ----------------------------------------------------------------------

def max_matching_convex(adj_H: list[list[int]], V: int):
    """
    Максимальне пароспівставлення в конвексному дводольному графі.

    Ліва доля: H = {0, ..., H-1}.
    Права доля: V = {0, ..., V-1}.

    Припущення: для кожного h сусіди adj_H[h] утворюють інтервал
    по індексах вертикалей (після сортування).

    Алгоритм:
      1. Для кожного h обчислюємо інтервал [L_h, R_h].
      2. Розглядаємо вертикалі v = 0..V-1 як "time slots".
      3. Для кожного v запускаємо greedy-алгоритм інтервального
         розкладу (unit jobs with release & deadline) з пріоритетом
         за мінімальним R_h.

    Час: O((H + V) log H) для самого matching; побудова adj_H
    через sweep займає O(E log H), де E – число перетинів.
    """
    H = len(adj_H)

    intervals = []
    for h_id, nbrs in enumerate(adj_H):
        if not nbrs:
            continue
        uniq = sorted(set(nbrs))
        L = uniq[0]
        R = uniq[-1]
        # Якщо хочете перевірити конвексність – розкоментуйте:
        # if R - L + 1 != len(uniq):
        #     raise ValueError(f"Non-convex adjacency for h={h_id}: {uniq}")
        intervals.append({"h_id": h_id, "L": L, "R": R})

    # Сортуємо інтервали за лівими кінцями
    intervals.sort(key=lambda it: it["L"])

    mate_H = [-1] * H
    mate_V = [-1] * V

    # Черга пріоритетів за правими кінцями R
    pq: list[tuple[int, int]] = []   # (R, h_id)
    p = 0  # покажчик по intervals

    for v_id in range(V):
        # Додаємо всі горизонталі, які вже "активні" (L_h <= v_id)
        while p < len(intervals) and intervals[p]["L"] <= v_id:
            it = intervals[p]
            if it["R"] >= v_id:
                heapq.heappush(pq, (it["R"], it["h_id"]))
            p += 1

        # Викидаємо прострочені (R_h < v_id)
        while pq and pq[0][0] < v_id:
            heapq.heappop(pq)

        # Беремо ту горизонталь, у якої R_h мінімальне
        if pq:
            R_best, h_id = heapq.heappop(pq)
            if mate_H[h_id] == -1 and mate_V[v_id] == -1:
                mate_H[h_id] = v_id
                mate_V[v_id] = h_id

    return mate_H, mate_V

# ----------------------------------------------------------------------
# 1. From maximum matching to maximum independent set in bipartite graph
# ----------------------------------------------------------------------

def max_independent_set_from_matching(adj_H, adj_V, mate_H, mate_V):
    """
    Обчислює максимальну незалежну множину вершин дводольного графа
    G = (H ∪ V, E) на основі максимального пароспівставлення M.

    Вхід:
        adj_H[h] - список v, суміжних з h (ліва доля H).
        adj_V[v] - список h, суміжних з v (права доля V).
        mate_H[h] - або v, з яким h зведено в M, або -1.
        mate_V[v] - або h, з яким v зведено в M, або -1.

    Вихід:
        I_H[h] = True, якщо h входить у максимальну незалежну множину.
        I_V[v] = True, якщо v входить у максимальну незалежну множину.

    Стандартний алгоритм:
        1) Побудувати мінімальне вершинне покриття C через пошук
           чергуючих шляхів з усіх вільних вершин лівої долі.
        2) Незалежна множина = доповнення до C.
    """
    H = len(adj_H)
    V = len(adj_V)

    # Z_H, Z_V — множини вершин, досяжних по чергуючих шляхах
    Z_H = [False] * H
    Z_V = [False] * V

    q = deque()

    # 1. Стартуємо з усіх вільних вершин лівої частини
    for h in range(H):
        if mate_H[h] == -1:
            Z_H[h] = True
            q.append(("H", h))

    # 2. BFS по чергуючих шляхах:
    #    H -> V лише по НЕспарених ребрах,
    #    V -> H лише по спарених.
    while q:
        side, u = q.popleft()
        if side == "H":
            h = u
            for v in adj_H[h]:
                if mate_H[h] == v:
                    # Спарене ребро — його використовуємо лише
                    # при переході V -> H
                    continue
                if not Z_V[v]:
                    Z_V[v] = True
                    q.append(("V", v))
        else:
            v = u
            h = mate_V[v]
            if h != -1 and not Z_H[h]:
                Z_H[h] = True
                q.append(("H", h))

    # Мінімальне вершинне покриття:
    #   C_H = H \ Z_H
    #   C_V = Z_V
    C_H = [not Z_H[h] for h in range(H)]
    C_V = [Z_V[v] for v in range(V)]

    # Максимальна незалежна множина = доповнення до покриття
    I_H = [not C_H[h] for h in range(H)]
    I_V = [not C_V[v] for v in range(V)]

    return I_H, I_V


# ----------------------------------------------------------------------
# 2. Додавання сегментів для "невкритих" вігнутих вершин (крок 3 Оhtsuki)
# ----------------------------------------------------------------------

def augment_with_remaining_reflexes(reflex_vertices,
                                    all_H, all_V,
                                    chosen_H, chosen_V):
    """
    Після вибору незалежної множини хорд (chosen_H, chosen_V),
    для кожної вігнутої вершини перевіряємо, чи є хоч один сегмент,
    що від неї виходить. Якщо ні — додаємо будь-який доступний
    горизонтальний або вертикальний відрізок видимості.

    Це відповідає ідеї кроку 3 у Ohtsuki: добудувати максимальні
    вертикальні/горизонтальні лінії від решти вігнутих вершин,
    щоб завершити розбиття.
    """
    used = {rv["id"]: False for rv in reflex_vertices}

    for seg in chosen_H:
        used[seg["vertex_id"]] = True
    for seg in chosen_V:
        used[seg["vertex_id"]] = True

    # Індексування всіх сегментів за vertex_id
    H_by_v = {}
    V_by_v = {}
    for seg in all_H:
        H_by_v.setdefault(seg["vertex_id"], []).append(seg)
    for seg in all_V:
        V_by_v.setdefault(seg["vertex_id"], []).append(seg)

    for rv in reflex_vertices:
        vid = rv["id"]
        if not used.get(vid, False):
            # Намагаємося спочатку використати горизонтальний,
            # якщо немає — вертикальний
            candidates = H_by_v.get(vid)
            if not candidates:
                candidates = V_by_v.get(vid)
            if not candidates:
                # Теоретично не повинно статись, але на всяк випадок
                continue
            seg = candidates[0]
            if "y" in seg:  # горизонтальний
                chosen_H.append(seg)
            else:          # вертикальний
                chosen_V.append(seg)
            used[vid] = True

    return chosen_H, chosen_V


# ----------------------------------------------------------------------
# 3. Візуалізація розбиття
# ----------------------------------------------------------------------

def draw_dissection(outer, holes, horiz_segments, vert_segments):
    """
    Малює:
      - контур зовнішнього багатокутника,
      - контури отворів,
      - обрані горизонтальні та вертикальні відрізки (хорди).

    Формат вершин: [x, y, is_reflex_bool]
    """
    fig, ax = plt.subplots()

    # Зовнішній контур
    xs = [v[0] for v in outer] + [outer[0][0]]
    ys = [v[1] for v in outer] + [outer[0][1]]
    ax.plot(xs, ys, "k-")

    # Отвори
    for hole in holes:
        xs = [v[0] for v in hole] + [hole[0][0]]
        ys = [v[1] for v in hole] + [hole[0][1]]
        ax.plot(xs, ys, "k-")

    # Горизонтальні відрізки
    for seg in horiz_segments:
        x1, x2 = seg["x_from"], seg["x_to"]
        y = seg["y"]
        ax.plot([x1, x2], [y, y], linewidth=1)

    # Вертикальні відрізки
    for seg in vert_segments:
        x = seg["x"]
        y1, y2 = seg["y_from"], seg["y_to"]
        ax.plot([x, x], [y1, y2], linewidth=1)

    ax.set_aspect("equal", "box")
    #ax.invert_yaxis()  # якщо координати у форматі екранних
    ax.grid(False)
    plt.show()


# ----------------------------------------------------------------------
# 4. Повна функція: від багатокутника до розбиття
# ----------------------------------------------------------------------

def partition_and_draw(outer, holes):
    """
    Повний цикл:
      1) Класифікуємо ребра, знаходимо вігнуті вершини.
      2) Будуємо горизонтальні та вертикальні відрізки видимості
         (\"хорди\") за допомогою sweep.
      3) Будуємо дводольний граф перетинів.
      4) Знаходимо максимальне пароспівставлення у конвексному
         дводольному графі (Imai–Asano).
      5) З нього виводимо максимальну незалежну множину хорд.
      6) Добудовуємо сегменти для решти вігнутих вершин.
      7) Малюємо кінцеве розбиття.

    Повертає:
      chosen_H, chosen_V — списки горизонтальних та вертикальних сегментів,
      які разом із контуром дають розбиття.
    """
    # 1. Передобробка (припускаємо, що вже виконано fix_orientation,
    #    determine_convexity, і outer/holes у форматі [x, y, is_reflex_bool])
    vertical_edges, horizontal_edges, reflex_vertices = classify_edges(
        outer, holes
    )

    # 2. Sweep-и: будуємо хорди
    horiz_chords = horizontal_sweep(vertical_edges, reflex_vertices)
    vert_chords = vertical_sweep(horizontal_edges, reflex_vertices)

    # 3. Дводольний граф перетинів хорд
    adj_H, adj_V = build_bipartite_graph(horiz_chords, vert_chords)

    # 4. Максимальне пароспівставлення у конвексному графі
    mate_H, mate_V = max_matching_convex(adj_H, len(vert_chords))

    # 5. Максимальна незалежна множина хорд
    I_H, I_V = max_independent_set_from_matching(adj_H, adj_V, mate_H, mate_V)

    chosen_H = [seg for h_id, seg in enumerate(horiz_chords) if I_H[h_id]]
    chosen_V = [seg for v_id, seg in enumerate(vert_chords) if I_V[v_id]]

    # 6. Добудувати сегменти від решти вігнутих вершин
    chosen_H, chosen_V = augment_with_remaining_reflexes(
        reflex_vertices,
        horiz_chords,
        vert_chords,
        chosen_H,
        chosen_V
    )

    # 7. Візуалізація
    draw_dissection(outer, holes, chosen_H, chosen_V)

    return chosen_H, chosen_V



def compute_partition(outer, holes):
    """
    Run the full Imai–Asano-style pipeline on the current polygon.

    `outer` and `holes` are lists of vertices [x, y] (before convexity).
    This function may modify them in-place to [x, y, is_reflex].

    Returns:
        outer, holes, chosen_H, chosen_V
    """
    # Preprocess outer
    validate(outer)
    fix_orientation(outer, is_hole=False)
    determine_convexity(outer, is_hole=False)

    # Preprocess holes
    for hole in holes:
        validate(hole)
        fix_orientation(hole, is_hole=True)
        determine_convexity(hole, is_hole=True)

    # classify_edges expects vertices as [x, y, is_reflex]
    vertical_edges, horizontal_edges, reflex_vertices = classify_edges(
        outer, holes
    )

    # Build visibility chords
    horiz_chords = horizontal_sweep(vertical_edges, reflex_vertices)
    vert_chords = vertical_sweep(horizontal_edges, reflex_vertices)

    # Bipartite intersection graph
    adj_H, adj_V = build_bipartite_graph(horiz_chords, vert_chords)

    # Maximum matching in convex bipartite graph
    mate_H, mate_V = max_matching_convex(adj_H, len(vert_chords))

    # Maximum independent set of chords
    I_H, I_V = max_independent_set_from_matching(adj_H, adj_V, mate_H, mate_V)

    chosen_H = [seg for h_id, seg in enumerate(horiz_chords) if I_H[h_id]]
    chosen_V = [seg for v_id, seg in enumerate(vert_chords) if I_V[v_id]]

    # Augment with segments from remaining reflex vertices
    chosen_H, chosen_V = augment_with_remaining_reflexes(
        reflex_vertices,
        horiz_chords,
        vert_chords,
        chosen_H,
        chosen_V,
    )

    return outer, holes, chosen_H, chosen_V


class ImaiAsanoGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("Imai–Asano Rectilinear Polygon Partition")

        # Data
        self.outer = []          # list[[x, y], ...]
        self.holes = []          # list[list[[x, y], ...], ...]
        self.chosen_H = []       # horizontal partition segments
        self.chosen_V = []       # vertical partition segments

        # Layout: top frame for buttons, bottom frame for canvas
        control_frame = tk.Frame(master)
        control_frame.pack(side=tk.TOP, fill=tk.X)

        canvas_frame = tk.Frame(master)
        canvas_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        # Buttons
        btn_load = tk.Button(
            control_frame,
            text="Load polygon from file",
            command=self.load_polygon_from_file,
        )
        btn_load.pack(side=tk.LEFT, padx=5, pady=5)

        btn_random = tk.Button(
            control_frame,
            text="Random rectilinear polygon",
            command=self.generate_random_polygon,
        )
        btn_random.pack(side=tk.LEFT, padx=5, pady=5)

        btn_partition = tk.Button(
            control_frame,
            text="Partition",
            command=self.partition_current_polygon,
        )
        btn_partition.pack(side=tk.LEFT, padx=5, pady=5)

        btn_clear = tk.Button(
            control_frame,
            text="Clear",
            command=self.clear_all,
        )
        btn_clear.pack(side=tk.LEFT, padx=5, pady=5)

        # Matplotlib figure embedded in Tkinter
        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect("equal", "box")
        self.ax.grid(False)

        self.canvas = FigureCanvasTkAgg(self.fig, master=canvas_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)

        self.redraw()

    # --------------------------------------------------------------
    # Polygon input
    # --------------------------------------------------------------

    def load_polygon_from_file(self):
        filename = filedialog.askopenfilename(
            title="Select polygon file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if not filename:
            return

        try:
            with open(filename, "r") as f:
                contents = f.read()
        except OSError as e:
            messagebox.showerror("Error", f"Cannot open file:\n{e}")
            return

        if "outer: " not in contents:
            messagebox.showerror(
                "Error",
                'File must contain outer polygon as\nouter: [[x1, y1], ..., [xn, yn]]'
            )
            return

        outer = []
        holes = []
        scanning_holes = False
        holes_description = ""

        for line in contents.splitlines():
            if "outer: " in line:
                content = line.replace("outer: ", "")
                outer = ast.literal_eval(content)
            elif "holes: " in line:
                scanning_holes = True
                content = line.replace("holes: ", "")
                holes_description += content
            elif scanning_holes:
                holes_description += line

        if holes_description.strip():
            holes = ast.literal_eval(holes_description)
        else:
            holes = []

        # Ensure we have [x, y] lists
        outer = [list(v) for v in outer]
        holes = [[list(v) for v in hole] for hole in holes]

        self.outer = outer
        self.holes = holes
        self.chosen_H = []
        self.chosen_V = []

        self.redraw()

    def generate_random_polygon(self, n_cols: int = 6,
                                width: int = 10, height: int = 10):
        """
        Generate a simple x-monotone rectilinear polygon (histogram-like)
        without holes, on a grid [0,width] x [0,height].
        """
        if n_cols < 2:
            n_cols = 2

        # Random x coordinates (columns)
        xs_inner = sorted(random.sample(range(1, width), n_cols - 1))
        xs = [0] + xs_inner + [width]

        # Random heights for each column break
        hs = [random.randint(1, height) for _ in xs]

        # Build upper chain: (x0,h0)->(x1,h0)->(x1,h1)->(x2,h1)->...
        verts = []
        verts.append([xs[0], 0])  # start on bottom

        for i in range(len(xs)):
            verts.append([xs[i], hs[i]])  # up
            if i < len(xs) - 1:
                verts.append([xs[i + 1], hs[i]])  # right

        verts.append([xs[-1], 0])  # down to bottom
        # Optionally remove consecutive duplicates
        simplified = []
        for v in verts:
            if not simplified or simplified[-1] != v:
                simplified.append(v)

        self.outer = simplified
        self.holes = []
        self.chosen_H = []
        self.chosen_V = []
        self.redraw()

    # --------------------------------------------------------------
    # Partition and drawing
    # --------------------------------------------------------------

    def partition_current_polygon(self):
        if not self.outer:
            messagebox.showinfo("Info", "No polygon loaded or generated.")
            return

        try:
            outer = [v[:] for v in self.outer]          # shallow copy
            holes = [[p[:] for p in hole] for hole in self.holes]

            outer, holes, chosen_H, chosen_V = compute_partition(outer, holes)

        except Exception as e:
            messagebox.showerror("Error during partition", str(e))
            return

        # Store results
        self.outer = outer
        self.holes = holes
        self.chosen_H = chosen_H
        self.chosen_V = chosen_V

        self.redraw()

    def clear_all(self):
        self.outer = []
        self.holes = []
        self.chosen_H = []
        self.chosen_V = []
        self.redraw()

    def redraw(self):
        self.ax.clear()
        self.ax.set_aspect("equal", "box")
        self.ax.grid(False)

        # Draw outer polygon
        if self.outer:
            xs = [v[0] for v in self.outer] + [self.outer[0][0]]
            ys = [v[1] for v in self.outer] + [self.outer[0][1]]
            self.ax.plot(xs, ys, "k-")

        # Draw holes
        for hole in self.holes:
            xs = [v[0] for v in hole] + [hole[0][0]]
            ys = [v[1] for v in hole] + [hole[0][1]]
            self.ax.plot(xs, ys, "k-")

        # Draw horizontal chords
        for seg in self.chosen_H:
            x1, x2 = seg["x_from"], seg["x_to"]
            y = seg["y"]
            self.ax.plot([x1, x2], [y, y], linewidth=1)

        # Draw vertical chords
        for seg in self.chosen_V:
            x = seg["x"]
            y1, y2 = seg["y_from"], seg["y_to"]
            self.ax.plot([x, x], [y1, y2], linewidth=1)

        self.canvas.draw()


if __name__ == "__main__":
    root = tk.Tk()
    app = ImaiAsanoGUI(root)
    root.mainloop()


if __name__ == '__main__':
    outer = []
    holes = []

    with open("input.txt", "r") as file:
        scanning_holes = False
        contents = file.read()

        holes_description = ""

        if "outer: " not in contents:
            raise InvalidPolygonError(
                'Файл мусить містити зовнішню оболонку "outer: [список вершин]"'
            )

        for line in contents.splitlines():
            if "outer: " in line:
                content = line.replace("outer: ", "")
                outer = ast.literal_eval(content)
            elif "holes: " in line:
                scanning_holes = True
                content = line.replace("holes: ", "")
                holes_description += content
            elif scanning_holes:
                holes_description += line

        if holes_description.strip():
            holes = ast.literal_eval(holes_description)
        else:
            holes = []

    # --- Preprocessing: validate, orient, convexity/reflexity ---

    validate(outer)
    fix_orientation(outer, is_hole=False)
    determine_convexity(outer, is_hole=False)

    for i, hole in enumerate(holes):
        validate(hole)
        fix_orientation(hole, is_hole=True)
        determine_convexity(hole, is_hole=True)

    # --- Run the partition algorithm and draw the result ---

    chosen_H, chosen_V = partition_and_draw(outer, holes)

