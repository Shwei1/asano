import math
from preprocess import *
from sortedcontainers import SortedKeyList, SortedList
import heapq
from collections import deque
from bisect import bisect_left


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
        events.append((e["x_high"], EDGE_END,   e["id"]))
        events.append((e["x_low"],  EDGE_START, e["id"]))

    for v in reflex_vertices:
        events.append((v["x"], REFLEX, v["id"]))

    events.sort(key=lambda ev: (ev[0], ev[1]))
    return events



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

    active_edges = SortedKeyList(key=lambda e: e["x"])
    
   
def horizontal_sweep(vertical_edges: list, reflex_vertices: list) -> list[dict]:
    events = build_events_y(vertical_edges, reflex_vertices)

    active_edges = SortedKeyList(key=lambda e: e["x"])
    edge_by_id = {e["id"]: e for e in vertical_edges}
    vtx_by_id = {v["id"]: v for v in reflex_vertices}

    visibility_segments: list[dict] = []

    for y, kind, obj_id in events:
        if kind == EDGE_START:
            active_edges.add(edge_by_id[obj_id])

        elif kind == EDGE_END:
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

    events = []

    for h_id, h in enumerate(horizontal_segs):
        events.append((h["_x_low"], 0, h_id))
        events.append((h["_x_high"], 2, h_id))

    vertical_segs_sorted = sorted(
        enumerate(vertical_segs),
        key=lambda kv: kv[1]["x"]
    )
    id_map = {}
    for new_id, (old_id, seg) in enumerate(vertical_segs_sorted):
        id_map[old_id] = new_id

    for old_id, seg in vertical_segs_sorted:
        x = seg["x"]
        v_id = id_map[old_id]
        events.append((x, 1, v_id))

    events.sort(key=lambda ev: (ev[0], ev[1]))

    active = SortedList()

    adj_H = [[] for _ in range(H)]
    adj_V = [[] for _ in range(V)]

    for x, kind, obj_id in events:
        if kind == 0:
            h_id = obj_id
            y = horizontal_segs[h_id]["y"]
            active.add((y, h_id))

        elif kind == 2:
            h_id = obj_id
            y = horizontal_segs[h_id]["y"]
            active.remove((y, h_id))

        else:
            v_id = obj_id
            if not active:
                continue

            old_id, seg = vertical_segs_sorted[v_id]
            v = seg
            y_low = v["_y_low"]
            y_high = v["_y_high"]

            lo = active.bisect_left((y_low, -math.inf))
            hi = active.bisect_right((y_high, math.inf))

            for idx in range(lo, hi):
                y, h_id = active[idx]
                adj_H[h_id].append(v_id)
                adj_V[v_id].append(h_id)

    return adj_H, adj_V


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
        intervals.append({"h_id": h_id, "L": L, "R": R})

    intervals.sort(key=lambda it: it["L"])

    mate_H = [-1] * H
    mate_V = [-1] * V

    pq: list[tuple[int, int]] = []  
    p = 0 

    for v_id in range(V):
        while p < len(intervals) and intervals[p]["L"] <= v_id:
            it = intervals[p]
            if it["R"] >= v_id:
                heapq.heappush(pq, (it["R"], it["h_id"]))
            p += 1

        while pq and pq[0][0] < v_id:
            heapq.heappop(pq)

        if pq:
            R_best, h_id = heapq.heappop(pq)
            if mate_H[h_id] == -1 and mate_V[v_id] == -1:
                mate_H[h_id] = v_id
                mate_V[v_id] = h_id

    return mate_H, mate_V


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

    Z_H = [False] * H
    Z_V = [False] * V

    q = deque()

    for h in range(H):
        if mate_H[h] == -1:
            Z_H[h] = True
            q.append(("H", h))

    while q:
        side, u = q.popleft()
        if side == "H":
            h = u
            for v in adj_H[h]:
                if mate_H[h] == v:
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

    C_H = [not Z_H[h] for h in range(H)]
    C_V = [Z_V[v] for v in range(V)]

    I_H = [not C_H[h] for h in range(H)]
    I_V = [not C_V[v] for v in range(V)]

    return I_H, I_V


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
            candidates = H_by_v.get(vid)
            if not candidates:
                candidates = V_by_v.get(vid)
            if not candidates:
                continue
            seg = candidates[0]
            if "y" in seg:  # горизонтальний
                chosen_H.append(seg)
            else:          # вертикальний
                chosen_V.append(seg)
            used[vid] = True

    return chosen_H, chosen_V

def compute_partition(outer, holes):
    """
    Run the full Imai–Asano-style pipeline on the current polygon.

    `outer` and `holes` are lists of vertices [x, y] (before convexity).
    This function may modify them in-place to [x, y, is_reflex].

    Returns:
        outer, holes, chosen_H, chosen_V
    """
    validate(outer)
    fix_orientation(outer, is_hole=False)
    determine_convexity(outer, is_hole=False)

    for hole in holes:
        validate(hole)
        fix_orientation(hole, is_hole=True)
        determine_convexity(hole, is_hole=True)

    vertical_edges, horizontal_edges, reflex_vertices = classify_edges(
        outer, holes
    )

    horiz_chords = horizontal_sweep(vertical_edges, reflex_vertices)
    vert_chords = vertical_sweep(horizontal_edges, reflex_vertices)

    adj_H, adj_V = build_bipartite_graph(horiz_chords, vert_chords)

    mate_H, mate_V = max_matching_convex(adj_H, len(vert_chords))

    I_H, I_V = max_independent_set_from_matching(adj_H, adj_V, mate_H, mate_V)

    chosen_H = [seg for h_id, seg in enumerate(horiz_chords) if I_H[h_id]]
    chosen_V = [seg for v_id, seg in enumerate(vert_chords) if I_V[v_id]]

    chosen_H, chosen_V = augment_with_remaining_reflexes(
        reflex_vertices,
        horiz_chords,
        vert_chords,
        chosen_H,
        chosen_V,
    )

    return outer, holes, chosen_H, chosen_V




