from exceptions import InvalidPolygonError

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
    Визначаємо опуклість вершин
    """
    n = len(v)
    for i in range(n):
        i_prev = (i - 1) % n
        i_next = (i + 1) % n

        vx_prev = v[i_prev]
        vx_curr = v[i]
        vx_next = v[i_next]

        x_prev, y_prev = vx_prev[0], vx_prev[1]
        x_cur,  y_cur  = vx_curr[0], vx_curr[1]
        x_next, y_next = vx_next[0], vx_next[1]

        prev_curr = (x_cur - x_prev, y_cur - y_prev)
        curr_next = (x_next - x_cur, y_next - y_cur)

        turn = prev_curr[0] * curr_next[1] - prev_curr[1] * curr_next[0]

        if is_hole:
            is_reflex = (turn < 0)
        else:
            is_reflex = (turn < 0)

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

            if x == x_next:  
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
            else:            
                x_low = min(x, x_next)
                x_high = max(x, x_next)
                left_to_right = (x_next > x)
                edge = {
                    "y": y, 
                    "x_low": x_low,
                    "x_high": x_high,
                    "loop_index": loop_index,
                    "is_outer": is_outer,
                    "i": i,
                    "left_to_right": left_to_right,
                }
                edge["id"] = len(horizontal_edges)
                horizontal_edges.append(edge)

    process_loop(outer, loop_index=0, is_outer=True)

    for h_idx, hole in enumerate(holes, start=1):
        process_loop(hole, loop_index=h_idx, is_outer=False)

    return vertical_edges, horizontal_edges, reflex_vertices

