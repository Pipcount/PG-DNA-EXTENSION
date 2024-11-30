from os import write

from pygments.styles.gh_dark import GRAY_3

from manimlib import *


class DataTypes(InteractiveScene):
    def construct(self) -> None:
        text_color = BLACK
        title_text_color = WHITE
        font = "Monospace"
        title_text_size = 100
        header_rect_color = "#032a3f"
        data_rectangle = Rectangle(width=FRAME_WIDTH, height=FRAME_HEIGHT, color=header_rect_color, fill_color=header_rect_color)
        data_rectangle.set_fill(header_rect_color, 1)
        self.add(data_rectangle)
        encoding_text = Text("Encoding", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        encoding_text.set_color(title_text_color)
        header_text = Text("Data Types", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        header_text.set_color(title_text_color)
        dna_header_text = Text("DNA", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        dna_header_text.set_color(title_text_color)
        kmer_header_text = Text("KMER", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        kmer_header_text.set_color(title_text_color)
        qkmer_header_text = Text("QKMER", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        qkmer_header_text.set_color(title_text_color)
        self.play(Write(header_text))
        self.wait(1)
        dna_header_text.to_edge(UP)
        kmer_header_text.to_edge(UP)
        qkmer_header_text.to_edge(UP)
        encoding_text.to_edge(UP)

        self.play(
            Transform(header_text, encoding_text),
            data_rectangle.animate.set_height(2.3, stretch=True).to_edge(UP).shift(UP),
            run_time=1
        )


        a = Text("A", font=font, font_size=80)
        c = Text("C", font=font, font_size=80)
        g = Text("G", font=font, font_size=80)
        t = Text("T", font=font, font_size=80)

        a.to_edge(LEFT)
        c.next_to(a, RIGHT, buff=0.5)  # c à droite de a
        g.next_to(c, RIGHT, buff=0.5)  # g à droite de c
        t.next_to(g, RIGHT, buff=0.5)  # t à droite de g

        translation_group = VGroup(a, c, g, t)
        translation_group.set_color(text_color)
        translation_group.move_to(ORIGIN)  # Déplace tout le groupe au centre de la scène

        self.play(Write(translation_group))
        self.wait(1)

        self.play(
            a.animate.shift(LEFT * 3),  # Déplacer A à gauche
            c.animate.shift(LEFT),  # Déplacer C à gauche
            g.animate.shift(RIGHT),  # Déplacer G à droite
            t.animate.shift(RIGHT * 3),  # Déplacer T à droite
            run_time=2  # Durée de l'animation
        )

        arrow_a = Tex(r"\rightarrow", font_size=72)
        arrow_c = Tex(r"\rightarrow", font_size=72)
        arrow_g = Tex(r"\rightarrow", font_size=72)
        arrow_t = Tex(r"\rightarrow", font_size=72)
        group_arrows = VGroup(arrow_a, arrow_c, arrow_g, arrow_t)
        group_arrows.set_color(text_color)

        arrow_a.next_to(a, RIGHT, buff=0.2)
        arrow_c.next_to(c, RIGHT, buff=0.2)
        arrow_g.next_to(g, RIGHT, buff=0.2)
        arrow_t.next_to(t, RIGHT, buff=0.2)

        # Ajouter les flèches à la scène
        self.play(Write(arrow_a), Write(arrow_c), Write(arrow_g), Write(arrow_t))

        # Créer les valeurs binaires
        binary_a = Tex("00", font_size=80)
        binary_c = Tex("01", font_size=80)
        binary_g = Tex("10", font_size=80)
        binary_t = Tex("11", font_size=80)
        group_binaries = VGroup(binary_a, binary_c, binary_g, binary_t)
        group_binaries.set_color(text_color)

        # Positionner les valeurs binaires après les flèches
        binary_a.next_to(arrow_a, RIGHT, buff=0.2)
        binary_c.next_to(arrow_c, RIGHT, buff=0.2)
        binary_g.next_to(arrow_g, RIGHT, buff=0.2)
        binary_t.next_to(arrow_t, RIGHT, buff=0.2)

        group_a = VGroup(a, arrow_a, binary_a)
        group_c = VGroup(c, arrow_c, binary_c)
        group_g = VGroup(g, arrow_g, binary_g)
        group_t = VGroup(t, arrow_t, binary_t)

        group_shift = VGroup(group_a, group_c, group_g, group_t)
        group_shift.generate_target()

        group_shift.target.move_to(ORIGIN)
        # Ajouter les valeurs binaires à la scène avec animation
        self.play(Write(binary_a), Write(binary_c), Write(binary_g), Write(binary_t),
                  MoveToTarget(group_shift), run_time=1)
        self.wait(1)

        self.play(
            group_a.animate.scale(0.5).to_edge(UP + LEFT).set_color(title_text_color),  # Déplacer le groupe A
            group_c.animate.scale(0.5).to_edge(UP + LEFT).shift(RIGHT * 1.3).set_color(title_text_color),  # Déplacer le groupe C
            group_g.animate.scale(0.5).to_edge(UP + LEFT).shift(DOWN * 0.5).set_color(title_text_color),  # Déplacer le groupe G
            group_t.animate.scale(0.5).to_edge(UP + LEFT).shift(DOWN * 0.5 + (RIGHT * 1.3)).set_color(title_text_color),
            run_time=1.5  # Durée de l'animation
        )

        self.wait(1)

        self.play(TransformMatchingStrings(header_text, dna_header_text, matched_keys=["D", "N", "A"]), run_time=1)

        dna_impl_text = Text("typedef bytea DNA;", font=font, font_size=50, t2w={"weight": BOLD}, )
        dna_impl_text.set_color(text_color)
        dna_impl_text.set_color_by_text("bytea", BLUE_E)
        dna_impl_text.move_to(dna_header_text.get_center() + 2 * DOWN)
        self.play(Write(dna_impl_text))

        self.wait(1)

        header_rectangle = Rectangle(width=(FRAME_WIDTH - 1) * 0.4, height=1.5, color=BLACK)
        data_rectangle = Rectangle(width=(FRAME_WIDTH - 1) * 0.6, height=1.5, color=BLACK)

        data_rectangle.next_to(header_rectangle, RIGHT, buff=0.01)

        rectangle_group = VGroup(header_rectangle, data_rectangle)
        rectangle_group.move_to(ORIGIN + 2 * DOWN)

        header_rectangle_text = Text("Header (4 Bytes)", font_size=50, color=WHITE)
        data_rectangle_text = Text("Data (n Bytes)", font_size=50, color=WHITE)
        header_rectangle_text.set_color(text_color)
        data_rectangle_text.set_color(text_color)

        header_rectangle_text.move_to(header_rectangle.get_center())
        data_rectangle_text.move_to(data_rectangle.get_center())

        self.play(Write(header_rectangle), Write(data_rectangle))
        self.play(Write(header_rectangle_text), Write(data_rectangle_text))
        data_rectangle_group = VGroup(data_rectangle_text, data_rectangle)
        data_rectangle_group.generate_target()
        data_rectangle_group.target.move_to(ORIGIN + 2 * DOWN)

        header_group = VGroup(header_rectangle_text, header_rectangle)

        self.play(
            FadeOut(header_group),
            FadeOut(dna_impl_text), MoveToTarget(data_rectangle_group)
        )
        self.play(data_rectangle.animate.set_width(FRAME_WIDTH - 1, stretch=True))

        dna_example_sequence = Text("ACGTCA", font=font, font_size=80)
        dna_example_sequence.set_color(text_color)
        dna_example_sequence.next_to(data_rectangle_text, 4 * UP, buff=0.5)
        self.play(Write(dna_example_sequence))

        translated_dna_example_sequence = Text("000110110100", font=font, font_size=70)
        translated_dna_example_sequence.set_color(text_color)
        translated_dna_example_sequence.next_to(dna_example_sequence, DOWN, buff=0.5)

        for i in range(6):
            self.play(
                dna_example_sequence[i].animate.set_color(PURPLE_E),
                TransformFromCopy(dna_example_sequence[i], translated_dna_example_sequence[2 * i:2 * i + 2]),
                run_time=1
            )
            self.wait(0.5)

        self.play(FadeOut(data_rectangle_text))

        last_byte_length_rect = Rectangle(width=data_rectangle.get_width()/3, height=data_rectangle.get_height(), color=data_rectangle.get_color())
        first_byte_rect = Rectangle(width=data_rectangle.get_width()/3, height=data_rectangle.get_height(), color=data_rectangle.get_color())
        second_byte_rect = Rectangle(width=data_rectangle.get_width()/3, height=data_rectangle.get_height(), color=data_rectangle.get_color())

        last_byte_length_rect.move_to(data_rectangle.get_center() + LEFT * data_rectangle.get_width()/3)
        first_byte_rect.move_to(data_rectangle.get_center())
        second_byte_rect.move_to(data_rectangle.get_center() + RIGHT * data_rectangle.get_width()/3)

        self.play(Write(last_byte_length_rect), Write(first_byte_rect), Write(second_byte_rect))

        first_byte_content = Text("00011011", font=font, font_size=70)
        first_byte_content.set_color(text_color)
        first_byte_content.move_to(first_byte_rect.get_center())

        second_byte_content = Text("0100", font=font, font_size=70)
        padding = Text("0000", font=font, font_size=70, t2c={"0": GRAY_3})
        padding.set_stroke(GRAY_3, 1)
        second_byte_content.set_color(text_color)
        second_byte_group = VGroup(second_byte_content, padding)
        second_byte_group.arrange(RIGHT, buff=0.2)
        second_byte_group.move_to(second_byte_rect.get_center())

        last_byte_length_content_0 = Text("00000000", font=font, font_size=70)
        last_byte_length_content_0.set_color(text_color)
        last_byte_length_content_0.move_to(last_byte_length_rect.get_center())

        self.play(TransformFromCopy(translated_dna_example_sequence[0:8], first_byte_content))
        self.play(TransformFromCopy(translated_dna_example_sequence[8:12], second_byte_content))
        self.play(Write(padding))
        self.play(TransformFromCopy(translated_dna_example_sequence[12:], last_byte_length_content_0))

        scanning_rect = SurroundingRectangle(second_byte_content[0:2], buff=0.1)
        scanning_rect.set_color(header_rect_color)

        self.play(ShowCreation(scanning_rect))
        updated_last_char = Text("1", font=font, font_size=70, color=YELLOW)
        updated_last_char.set_color(BLACK)
        updated_last_char.move_to(last_byte_length_content_0[-1].get_center())  # Position the new last character

        # Animate the transformation of the last character
        self.play(Transform(last_byte_length_content_0[-1], updated_last_char))
        self.play(scanning_rect.animate.surround(second_byte_content[2:4]))

        # Second animation: Transform second last '0' to '1' and last '1' to '0'
        second_last_char = last_byte_length_content_0[-2]
        updated_second_last_char = Text("1", font=font, font_size=70)
        updated_second_last_char.set_color(BLACK)
        updated_second_last_char.move_to(second_last_char.get_center())

        last_char = last_byte_length_content_0[-1]
        updated_last_char = Text("0", font=font, font_size=70)
        updated_last_char.set_color(BLACK)
        updated_last_char.move_to(last_char.get_center())

        # Animate the transformations of both the last and second-to-last characters
        self.play(
            Transform(second_last_char, updated_second_last_char),
            Transform(last_char, updated_last_char)
        )
        self.play(FadeOut(scanning_rect))
