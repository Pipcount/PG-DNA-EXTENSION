from manimlib import *

class test(InteractiveScene):
    def construct(self) -> None:
        text_color = BLACK
        font = "Monospace"
        encoding_text = Text("Encoding", font=font, font_size=80)
        encoding_text.set_color(text_color)
        header_text = Text("Data Types", font_size=80)
        header_text.set_color(text_color)
        dna_header_text = Text("DNA", font=font, font_size=80)
        dna_header_text.set_color(text_color)
        kmer_header_text = Text("KMER", font=font, font_size=80)
        kmer_header_text.set_color(text_color)
        qkmer_header_text = Text("QKMER", font=font, font_size=80)
        qkmer_header_text.set_color(text_color)
        self.play(Write(header_text))
        self.wait(1)
        dna_header_text.to_edge(UP)
        kmer_header_text.to_edge(UP)
        qkmer_header_text.to_edge(UP)
        encoding_text.to_edge(UP)

        self.play(Transform(header_text, encoding_text), run_time=1)

        a = Text("A", font_size=80, color=BLUE)
        c = Text("C", font_size=80, color=RED)
        g = Text("G", font_size=80, color=GREEN)
        t = Text("T", font_size=80, color=YELLOW)

        a.to_edge(LEFT)
        c.next_to(a, RIGHT, buff=0.5)  # c à droite de a
        g.next_to(c, RIGHT, buff=0.5)  # g à droite de c
        t.next_to(g, RIGHT, buff=0.5)  # t à droite de g

        group = VGroup(a, c, g, t)
        group.set_color(text_color)
        group.move_to(ORIGIN)  # Déplace tout le groupe au centre de la scène

        self.play(Write(group))
        self.wait(1)

        self.play(
            a.animate.shift(LEFT * 3),    # Déplacer A à gauche
            c.animate.shift(LEFT ),        # Déplacer C à gauche
            g.animate.shift(RIGHT),       # Déplacer G à droite
            t.animate.shift(RIGHT * 3),   # Déplacer T à droite
            run_time=2                    # Durée de l'animation
        )

        arrow_a = Text("→", font_size=72, color=WHITE)
        arrow_c = Text("→", font_size=72, color=WHITE)
        arrow_g = Text("→", font_size=72, color=WHITE)
        arrow_t = Text("→", font_size=72, color=WHITE)
        group_arrows = VGroup(arrow_a, arrow_c, arrow_g, arrow_t)
        group_arrows.set_color(text_color)

        arrow_a.next_to(a, RIGHT, buff=0.2)
        arrow_c.next_to(c, RIGHT, buff=0.2)
        arrow_g.next_to(g, RIGHT, buff=0.2)
        arrow_t.next_to(t, RIGHT, buff=0.2)

        # Ajouter les flèches à la scène
        self.play(Write(arrow_a), Write(arrow_c), Write(arrow_g), Write(arrow_t))

        # Créer les valeurs binaires
        binary_a = Text("00", font_size=80, color=BLUE_E)
        binary_c = Text("01", font_size=80, color=RED_E)
        binary_g = Text("10", font_size=80, color=GREEN_E)
        binary_t = Text("11", font_size=80, color=YELLOW_E)
        group_binaries = VGroup(binary_a, binary_c, binary_g, binary_t)
        group_binaries.set_color(text_color)

        # Positionner les valeurs binaires après les flèches
        binary_a.next_to(arrow_a, RIGHT, buff=0.2)
        binary_c.next_to(arrow_c, RIGHT, buff=0.2)
        binary_g.next_to(arrow_g, RIGHT, buff=0.2)
        binary_t.next_to(arrow_t, RIGHT, buff=0.2)

        # Ajouter les valeurs binaires à la scène avec animation
        self.play(Write(binary_a), Write(binary_c), Write(binary_g), Write(binary_t))
        self.wait(2)

        group_a = VGroup(a, arrow_a, binary_a)
        group_c = VGroup(c, arrow_c, binary_c)
        group_g = VGroup(g, arrow_g, binary_g)
        group_t = VGroup(t, arrow_t, binary_t)

        self.play(
            group_a.animate.scale(0.5).to_edge(UP + LEFT),  # Déplacer le groupe A
            group_c.animate.scale(0.5).to_edge(UP + LEFT).shift(RIGHT * 1.3),  # Déplacer le groupe C
            group_g.animate.scale(0.5).to_edge(UP + LEFT).shift(DOWN * 0.5),  # Déplacer le groupe G
            group_t.animate.scale(0.5).to_edge(UP + LEFT).shift(DOWN * 0.5 + (RIGHT * 1.3)), 
            run_time=3  # Durée de l'animation
        )

        # self.wait(1)
        group_1 = VGroup(group_a, group_c)
        group_2 = VGroup(group_g, group_t)
        total_group = VGroup(group_1, group_2)

        # Créer un grand rectangle qui englobe tout le groupe
        enclosing_rect = Rectangle(
            color=BLACK,
            width=total_group.get_width() + 0.7,  # Ajouter un petit espacement autour du groupe
            height=total_group.get_height() + 0.5  # Ajouter un petit espacement autour du groupe
        )

        # Centrer le rectangle autour du groupe
        enclosing_rect.move_to(total_group.get_center())

        # Ajouter le rectangle à la scène
        self.play(
            Write(enclosing_rect) 
        )

        
        self.wait(1)

        self.play(TransformMatchingStrings(header_text, dna_header_text, matched_keys=["D", "N", "A"]), run_time=1)

        dna_impl_text = Text("typedef bytea DNA;", font=font, font_size=50, t2w={"weight": BOLD},)
        dna_impl_text.set_color(text_color)
        dna_impl_text.set_color_by_text("bytea", BLUE_E)
        dna_impl_text.move_to(dna_header_text.get_center() + 2 * DOWN)
        self.play(Write(dna_impl_text))            
        
        self.wait(1)