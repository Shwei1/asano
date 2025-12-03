from exceptions import InvalidPolygonError
from sortedcontainers import SortedList
import ast

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

    

def determine_convexity(v: list, is_hole=False) -> None:
    """
    Визначення опуклості або вігнутості багатокутника в усіх вершинах
    Якщо вершині приписується True, то вона вігнута
    Якщо True, то опукла.
    """
    n = len(v)
    for i in range(n):
        i_prev = (i - 1) % n
        i_next = (i + 1) % n
        vertex_prev = v[i_prev]
        vertex_curr = v[i]
        vertex_next = v[i_next]

        prev_curr = [vertex_curr[0] - vertex_prev[0], vertex_curr[1] - vertex_prev[1]]
        curr_next = [vertex_next[0] - vertex_curr[0], vertex_next[1] - vertex_curr[1]]

        turn = prev_curr[0]*curr_next[1]-prev_curr[1]*curr_next[0]

        if is_hole:
            if turn < 0:
                vertex_curr.append(False)
            else:
                vertex_curr.append(True)
        else:
            if turn > 0:
                vertex_curr.append(False)
            else:
                vertex_curr.append(True)


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


def build_events_y(edges: list[dict], vertices: list[dict]):
    """
    Визначення типів подій, що відбуваються вздовж певної осі
    0 - початок проходження ребра
    1 - знаходження вігнутої вершини 
    2 - кінець проходження ребра
    """
    events = []

    for edge in edges:
        events.append((edge["y_low"], 0, edge["id"]))
        events.append((edge["y_high"], 2, edge["id"]))

    for vertex in vertices:
        events.append((vertex["y"], 1, vertex["id"]))

    events.sort(key=lambda e: (e[0], e[1]))

        

def horizontal_sweep(vertical_edges: list, reflex_vertices: list) -> list:
    """
    Будує вертикальні відрізки видимості за допомогою горизонтального сканування
    """
    tree = SortedList()

    events = build_events(vertical_edges, reflex_vertices)

    edge_by_id = {e["id"]: e for e in vertical_edges}
    vtx_by_id = {v["id"]: v for v in reflex_vertices}

    visibility_segments = []

    for y, kind, obj_id in events:
        if kind == 0:
            e = edge_by_id[obj_id]
            tree.add(e)

        elif kind == 2:
            e = edge_by_id[obj_id]
            tree.remove(e)

        else:
            v = vtx_by_id[obj_id]
            x_v = v["x"]

            if len(tree) == 0:
                continue

            xs = [e["x"] for e in tree]
            m = len(xs)

            inside = [(i%2==0) for i in range(m-1)]

            j = tree.bisect_left(xs, x_v)

            if j < m and xs[j] == x_v:
                if j > 0 and inside[j-1]:
                    visibility_segments.append({
                        "y": y,
                        "x_from": x_v,
                        "x_to": xs[j - 1],
                        "vertex_id": v["id"],
                        "direction": "left",
                    })
                if j < m - 1 and inside[j]:
                    visibility_segments.append({
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
                            "y": y,
                            "x_from": x_v,
                            "x_to": xs[idx],     
                            "vertex_id": v["id"],
                            "direction": "left",
                        })
                        visibility_segments.append({
                            "y": y,
                            "x_from": x_v,
                            "x_to": xs[idx + 1],  
                            "vertex_id": v["id"],
                            "direction": "right",
                        })


    return visibility_segments


def vertical_sweep(horizontal_edges: list, reflex_vertices: list) -> list:
    """
    Аналог horizontal_sweep, але для вертикальних відрізків видимості.

    Використовує горизонтальні ребра як межі, сканує по x,
    і від кожної вігнутої вершини пускає вертикальні промені вгору/вниз,
    поки не дійде до найближчої горизонтальної межі.

    Повертає список вертикальних відрізків:
        segment = {
            "x": x,
            "y_from": y0,
            "y_to": y1,
            "vertex_id": v_id,
            "direction": "down" / "up",
        }
    """
    events = build_events_x(horizontal_edges, reflex_vertices)

    edge_by_id = {e["id"]: e for e in horizontal_edges}
    vtx_by_id = {v["id"]: v for v in reflex_vertices}

    # У дереві тримаємо y-координати активних горизонтальних ребер
    active = SortedList()

    visibility_segments = []

    for x, kind, obj_id in events:
        if kind == 0:
            # EDGE_START: ребро стає активним
            e = edge_by_id[obj_id]
            active.add(e["y"])

        elif kind == 2:
            # EDGE_END: ребро втрачає активність
            e = edge_by_id[obj_id]
            active.remove(e["y"])

        else:
            # REFLEX
            v = vtx_by_id[obj_id]
            y_v = v["y"]

            if not active:
                continue

            ys = list(active)
            m = len(ys)

            # inside[i] відповідає інтервалу (ys[i], ys[i+1])
            inside = [(i % 2 == 0) for i in range(m - 1)]

            j = active.bisect_left(y_v)

            if j < m and ys[j] == y_v:
                # Вершина лежить на горизонтальному ребрі y = ys[j]
                # Нижче: (ys[j-1], ys[j])
                if j > 0 and inside[j - 1]:
                    visibility_segments.append({
                        "x": x,
                        "y_from": y_v,
                        "y_to": ys[j - 1],
                        "vertex_id": v["id"],
                        "direction": "down",
                    })
                # Вище: (ys[j], ys[j+1])
                if j < m - 1 and inside[j]:
                    visibility_segments.append({
                        "x": x,
                        "y_from": y_v,
                        "y_to": ys[j + 1],
                        "vertex_id": v["id"],
                        "direction": "up",
                    })
            else:
                # Вершина між ys[j-1] та ys[j]
                if 0 < j < m:
                    idx = j - 1
                    if inside[idx]:
                        visibility_segments.append({
                            "x": x,
                            "y_from": y_v,
                            "y_to": ys[idx],      # вниз
                            "vertex_id": v["id"],
                            "direction": "down",
                        })
                        visibility_segments.append({
                            "x": x,
                            "y_from": y_v,
                            "y_to": ys[idx + 1],  # вгору
                            "vertex_id": v["id"],
                            "direction": "up",
                        })

    return visibility_segments


def build_bipartite_graph(horizontal_segs: list, vertical_segs: list):
    """
    Створює дводольний граф перетинів між горизонтальними та вертикальними
    відрізками видимості.

    Вхід:
        horizontal_segs: список словників { "y", "x_from", "x_to", ... }
        vertical_segs:   список словників { "x", "y_from", "y_to", ... }

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
        events.append((h["_x_low"],  0, h_id))
        events.append((h["_x_high"], 2, h_id))

    for v_id, v in enumerate(vertical_segs):
        events.append((v["x"], 1, v_id))

    # Порядок: H_START (0), V_QUERY (1), H_END (2)
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

            v = vertical_segs[v_id]
            y_low = v["_y_low"]
            y_high = v["_y_high"]

            # Знаходимо всі активні горизонтальні з y в [y_low, y_high]
            lo = active.bisect_left((y_low, -math.inf))
            hi = active.bisect_right((y_high, math.inf))

            for idx in range(lo, hi):
                y, h_id = active[idx]
                # На цій x ми вже гарантуємо, що горизонтальний сегмент
                # покриває x (бо він активний), і y лежить в [y_low, y_high],
                # отже є перетин.
                adj_H[h_id].append(v_id)
                adj_V[v_id].append(h_id)

    return adj_H, adj_V

   












if __name__ == '__main__':
    outer = []
    holes = []

    with open("input.txt", "r") as file:
        scanning_holes = False
        contents = file.read()

        holes_description = ""

        if "outer: " not in contents:
            raise InvalidPolygonError("Файл мусить містити зовнішню оболонку \"outer: [список вершин]\"")

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
        
        holes = ast.literal_eval(holes_description)

    validate(outer)
    fix_orientation(outer)
    determine_convexity(outer)

    for hole in holes:
        validate(hole)
        fix_orientation(hole)
        determine_convexity(hole)
        for vertex in hole:
            reflex = vertex[2]
            if reflex:
                print(vertex, "is reflex")




    for vertex in outer:
        reflex = vertex[2]
        if reflex:
            print(vertex, " is reflex")


                






            



