from sortedcontainers import SortedList, SortedKeyList

def draw_dissection(outer, holes, horiz_segments, vert_segments):
    """
    Малює:
      - контур зовнішнього багатокутника,
      - контури отворів,
      - обрані горизонтальні та вертикальні відрізки (хорди).

    Формат вершин: [x, y, is_reflex_bool]
    """
    fig, ax = plt.subplots()

    xs = [v[0] for v in outer] + [outer[0][0]]
    ys = [v[1] for v in outer] + [outer[0][1]]
    ax.plot(xs, ys, "k-")

    for hole in holes:
        xs = [v[0] for v in hole] + [hole[0][0]]
        ys = [v[1] for v in hole] + [hole[0][1]]
        ax.plot(xs, ys, "k-")

    for seg in horiz_segments:
        x1, x2 = seg["x_from"], seg["x_to"]
        y = seg["y"]
        ax.plot([x1, x2], [y, y], linewidth=1)

    for seg in vert_segments:
        x = seg["x"]
        y1, y2 = seg["y_from"], seg["y_to"]
        ax.plot([x, x], [y1, y2], linewidth=1)

    ax.set_aspect("equal", "box")
    ax.grid(False)
    plt.show()


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
        chosen_V
    )

