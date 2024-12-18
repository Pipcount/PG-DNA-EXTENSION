from idlelib.pyparse import trans
from os import write
from manimlib import *


class KmeAnim(InteractiveScene):
    def construct(self) -> None:
        text_color = BLACK
        title_text_color = WHITE
        header_rect_size = 1.9
        font = "Monospace"
        title_text_size = 100
        header_rect_color = "#032a3f"
        bullet_point_text_size = 40

        top_of_slide = TOP - header_rect_size * UP
        bottom_of_slide = BOTTOM
        middle_of_slide = (top_of_slide + bottom_of_slide) / 2
        bullet_point_pos = top_of_slide + DOWN + 3 * LEFT

        header_rectangle = Rectangle(width=FRAME_WIDTH, height=FRAME_HEIGHT, color=header_rect_color, fill_color=header_rect_color)
        header_rectangle.set_fill(header_rect_color, 1)
        self.add(header_rectangle)
        data_types_text = Text("Data Types", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        data_types_text.set_color(title_text_color)
        encoding_text = Text("Encoding", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        encoding_text.set_color(title_text_color)
        dna_header_text = Text("DNA", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        dna_header_text.set_color(title_text_color)
        kmer_header_text = Text("KMER", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        kmer_header_text.set_color(title_text_color)
        qkmer_header_text = Text("QKMER", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        qkmer_header_text.set_color(title_text_color)
        self.play(Write(data_types_text))
        self.wait(1)
        dna_header_text.to_edge(UP)
        kmer_header_text.to_edge(UP)
        qkmer_header_text.to_edge(UP)
        encoding_text.to_edge(UP)

        header_rect_pos = TOP - header_rect_size / 2 * UP
        self.play(
            TransformMatchingShapes(data_types_text, encoding_text),
            header_rectangle.animate.set_height(header_rect_size, stretch=True).move_to(header_rect_pos),
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

        self.play(TransformMatchingStrings(encoding_text, dna_header_text, matched_keys=["d", "n", "e"]), run_time=1)

        dna_impl_text = Text("typedef bytea DNA;", font=font, font_size=50, t2w={"weight": BOLD}, )
        dna_impl_text.set_color(text_color)
        dna_impl_text.set_color_by_text("bytea", BLUE_E)
        dna_impl_text.move_to(dna_header_text.get_center() + 2 * DOWN)
        self.play(Write(dna_impl_text))

        self.wait(1)

        dna_header_rectangle = Rectangle(width=(FRAME_WIDTH - 1) * 0.4, height=1.5, color=BLACK)
        data_rectangle = Rectangle(width=(FRAME_WIDTH - 1) * 0.6, height=1.5, color=BLACK)

        data_rectangle.next_to(dna_header_rectangle, RIGHT, buff=0.01)

        rectangle_group = VGroup(dna_header_rectangle, data_rectangle)
        rectangle_group.move_to(ORIGIN + 2 * DOWN)

        header_rectangle_text = Text("Header (4 Bytes)", font_size=50, color=WHITE)
        data_rectangle_text = Text("Data (n Bytes)", font_size=50, color=WHITE)
        header_rectangle_text.set_color(text_color)
        data_rectangle_text.set_color(text_color)

        header_rectangle_text.move_to(dna_header_rectangle.get_center())
        data_rectangle_text.move_to(data_rectangle.get_center())

        self.play(Write(dna_header_rectangle), Write(data_rectangle))
        self.play(Write(header_rectangle_text), Write(data_rectangle_text))
        data_rectangle_group = VGroup(data_rectangle_text, data_rectangle)
        data_rectangle_group.generate_target()
        data_rectangle_group.target.move_to(ORIGIN + 2 * DOWN)

        header_group = VGroup(header_rectangle_text, dna_header_rectangle)

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
        padding = Text("0000", font=font, font_size=70, t2c={"0": GREY_C})
        padding.set_stroke(GREY_C, 1)
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

        translated_example_group = VGroup(first_byte_content, second_byte_group, last_byte_length_content_0,
                                          first_byte_rect, second_byte_rect, last_byte_length_rect, data_rectangle)

        self.play(
            TransformMatchingStrings(dna_header_text, kmer_header_text),
            FadeOut(first_byte_content),
            FadeOut(second_byte_content),
            FadeOut(padding),
            FadeOut(last_byte_length_content_0),
            FadeOut(last_byte_length_rect),
            FadeOut(first_byte_rect),
            FadeOut(second_byte_rect),
            FadeOut(dna_example_sequence),
            FadeOut(translated_dna_example_sequence),
            FadeOut(data_rectangle),
            run_time=1
        )

        kmer_struct_text_wo_len = Text(
            """
            struct Kmer {
                uint64_t kmer;
            } Kmer;
            """, font=font, font_size=50, t2w={"weight": BOLD}
        )

        kmer_struct_text_wo_len.set_color(text_color)
        kmer_struct_text_wo_len.set_color_by_text("struct", BLUE_E)
        kmer_struct_text_wo_len.set_color_by_text("uint64_t", BLUE_E)
        kmer_struct_text_wo_len.set_color_by_text("Kmer", ORANGE)

        kmer_struct_text_with_len = Text(
            """
            struct Kmer {
                uint64_t kmer;
                uint8_t k;
            } Kmer;
            """, font=font, font_size=50, t2w={"weight": BOLD}
        )
        kmer_struct_text_with_len.set_color(text_color)
        kmer_struct_text_with_len.set_color_by_text("struct", BLUE_E)
        kmer_struct_text_with_len.set_color_by_text("uint64_t", BLUE_E)
        kmer_struct_text_with_len.set_color_by_text("Kmer", ORANGE)
        kmer_struct_text_with_len.set_color_by_text("uint8_t", BLUE_E)
        kmer_struct_text_with_len.move_to(middle_of_slide)

        kmer_struct_text_wo_len.align_to(kmer_struct_text_with_len, UP)

        self.play(Write(kmer_struct_text_wo_len))
        self.wait(1)

        self.play(TransformMatchingStrings(kmer_struct_text_wo_len, kmer_struct_text_with_len, matched_keys=["uint", "_t", "k", ";"]), run_time=1)
        self.wait(1)

        self.play(
            FadeOut(kmer_struct_text_with_len),
            TransformMatchingStrings(kmer_header_text, qkmer_header_text, matched_keys=["q", "k", "m", "e", "r"]),
            run_time=1
        )

        qkmer_struct_text = Text(
            """
            struct Qkmer {
                uint64_t ac;
                uint64_t gt;
                uint8_t k;
            } Qkmer;
            """, font=font, font_size=50, t2w={"weight": BOLD}
        )
        qkmer_struct_text.set_color(text_color)
        qkmer_struct_text.set_color_by_text("struct", BLUE_E)
        qkmer_struct_text.set_color_by_text("uint64_t", BLUE_E)
        qkmer_struct_text.set_color_by_text("Qkmer", ORANGE)
        qkmer_struct_text.set_color_by_text("uint8_t", BLUE_E)
        qkmer_struct_text.move_to(middle_of_slide)

        self.play(Write(qkmer_struct_text))
        self.wait(1)

        length_header_text = Text("Length functions", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        length_header_text.set_color(title_text_color)
        length_header_text.to_edge(UP)

        self.play(
            FadeOut(qkmer_struct_text),
            FadeOut(group_shift),
            TransformMatchingStrings(qkmer_header_text, length_header_text),
            run_time=1
        )

        dna_sequence_text = Text("DNA", font=font, font_size=50)
        dna_sequence_text.set_color(text_color)
        dna_sequence_text.set_stroke(text_color, 2)
        kmer_text = Text("Kmer", font=font, font_size=50)
        kmer_text.set_color(text_color)
        kmer_text.set_stroke(text_color, 2)
        qkmer_text = Text("Qkmer", font=font, font_size=50)
        qkmer_text.set_color(text_color)
        qkmer_text.set_stroke(text_color, 2)

        # types_group = VGroup(dna_sequence_text, kmer_text, qkmer_text)
        # types_group.arrange(DOWN, buff=1.5)
        # types_group.move_to(middle_of_slide + UP/2)

        dna_length_function_text = Text(
            """
            uint32_t length = (VARSIZE(dna) - VARHDRSZ - 1) * 4;
            return length - (4 - last_byte_length(dna));
            """, font=font, font_size=40, t2w={"weight": BOLD}
        )
        dna_length_function_text.set_color(text_color)
        dna_length_function_text.set_color_by_text("uint32_t", BLUE_E)
        dna_length_function_text.set_color_by_text("VARSIZE", ORANGE)
        dna_length_function_text.set_color_by_text("VARHDRSZ", ORANGE)
        dna_length_function_text.set_color_by_text("return", BLUE_E)
        dna_length_function_text.set_color_by_text("last_byte_length", BLUE_E)

        kmer_length_function_text = Text(
            """
            return kmer.k;
            """, font=font, font_size=40, t2w={"weight": BOLD}
        )
        kmer_length_function_text.set_color(text_color)
        kmer_length_function_text.set_color_by_text("return", BLUE_E)

        qkmer_length_function_text = Text(
            """
            return qkmer.k;
            """, font=font, font_size=40, t2w={"weight": BOLD}
        )
        qkmer_length_function_text.set_color(text_color)
        qkmer_length_function_text.set_color_by_text("return", BLUE_E)


        dna_group = VGroup(dna_sequence_text, dna_length_function_text)
        dna_group.arrange(DOWN, buff=0.2)

        kmer_group = VGroup(kmer_text, kmer_length_function_text)
        kmer_group.arrange(DOWN, buff=0.2)

        qkmer_group = VGroup(qkmer_text, qkmer_length_function_text)
        qkmer_group.arrange(DOWN, buff=0.2)

        length_functions_group = VGroup(dna_group, kmer_group, qkmer_group)
        length_functions_group.arrange(DOWN, buff=1)
        length_functions_group.move_to(middle_of_slide)

        self.play(
            Write(length_functions_group),
            run_time=1
        )

        functions_text = Text("Functions", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        functions_text.set_color(title_text_color)
        functions_text.move_to(ORIGIN)

        self.play(
            FadeOut(length_functions_group, run_time=0.5),
            TransformMatchingShapes(length_header_text, functions_text),
            header_rectangle.animate.set_height(FRAME_HEIGHT, stretch=True).move_to(ORIGIN),
            run_time=1
        )
        self.wait()

        generate_kmers_header = Text(
            """
            Generate Kmers
            """, font=font, font_size=title_text_size, t2w={"weight": BOLD}
        )
        generate_kmers_header.set_color(title_text_color)
        generate_kmers_header.to_edge(UP)

        self.play(
            TransformMatchingShapes(functions_text, generate_kmers_header),
            header_rectangle.animate.set_height(header_rect_size, stretch=True).move_to(header_rect_pos),
            run_time=1
        )

        # Choses à dire:
        # C'est une set returning function
        # Donc dans le code C, il faut faire la structure typique pour une fonction de ce type
        # On ne présente pas cette fonction
        generate_kmers_sql = Text(
            r"""
            CREATE OR REPLACE FUNCTION generate_kmers(DNA, integer)
            RETURNS SETOF kmer
            """, font=font, font_size=40, t2w={"weight": BOLD}
        )
        generate_kmers_sql.set_color(text_color)
        generate_kmers_sql.set_color_by_text("CREATE", PURPLE_E)
        generate_kmers_sql.set_color_by_text("OR", PURPLE_E)
        generate_kmers_sql.set_color_by_text("REPLACE", PURPLE_E)
        generate_kmers_sql.set_color_by_text("RETURNS", PURPLE_E)
        generate_kmers_sql.set_color_by_text("integer", BLUE_E)
        generate_kmers_sql.set_color_by_text("FUNCTION", ORANGE)
        generate_kmers_sql.set_color_by_text("SETOF", PURPLE_E)
        generate_kmers_sql.move_to(middle_of_slide)

        self.play(Write(generate_kmers_sql))

        self.wait()
        translated_example_group.move_to(middle_of_slide + UP)
        dna_example_sequence.next_to(translated_example_group, UP, buff=0.5)
        last_byte_length_rect.set_fill(ORANGE, 0.5)
        self.play(
            TransformMatchingShapes(generate_kmers_sql, dna_example_sequence),
            FadeIn(translated_example_group),
            run_time=2
        )

        byte_ptr = Text("byte pointer", font=font, font_size=40, t2w={"weight": BOLD})
        nucleotide_ctr_text = Text("nucleotide counter:", font=font, font_size=40, t2w={"weight": BOLD})
        nucleotide_ctr = Text("0", font=font, font_size=40, t2w={"weight": BOLD})
        kmer_length_info = Text("kmer length: 4", font=font, font_size=40, t2w={"weight": BOLD})
        kmer_length_info.set_color(text_color)

        output_kmer = Text("Output: ", font=font, font_size=60, t2w={"weight": BOLD})
        output_kmer.set_color(text_color)
        output_kmer.set_stroke(text_color, 2)

        byte_ptr_arrow = Arrow(DOWN, UP/3, color=text_color)

        byte_ptr_group = VGroup(byte_ptr, byte_ptr_arrow)
        nucleotide_ctr_group = VGroup(nucleotide_ctr_text, nucleotide_ctr)
        byte_ptr_group.arrange(UP, buff=0.2)
        nucleotide_ctr_group.arrange(RIGHT, buff=0.2)
        byte_ptr_group.set_color(text_color)
        info_group = VGroup(nucleotide_ctr_group, kmer_length_info)
        info_group.arrange(DOWN, buff=0.2)
        info_group.set_color(text_color)
        kmer_length_info.align_to(nucleotide_ctr_text, LEFT)

        byte_ptr_group.move_to(first_byte_content[0].get_center() + 3 * DOWN/ 2)
        info_group.to_corner(DL)
        output_kmer.next_to(info_group, RIGHT, buff=1.5)

        self.play(Write(byte_ptr_group), Write(info_group), Write(output_kmer))

        self.wait()

        output_kmer_i = Text("ACGT", font=font, font_size=60, t2w={"weight": BOLD})
        output_kmer_i.set_color(text_color)
        output_kmer_i.set_stroke(text_color, 2)
        output_kmer_i.next_to(output_kmer, RIGHT, buff=0.5)
        self.play(
            TransformFromCopy(first_byte_content[0:2], output_kmer_i[0]),
            run_time = 1
        )

        new_ctr = Text("1", font=font, font_size=40, t2w={"weight": BOLD})
        new_ctr.set_color(text_color)
        new_ctr.move_to(nucleotide_ctr.get_center())
        self.play(
            Transform(nucleotide_ctr, new_ctr),
            TransformFromCopy(first_byte_content[2:4], output_kmer_i[1]),
            run_time=1
        )

        new_ctr = Text("2", font=font, font_size=40, t2w={"weight": BOLD})
        new_ctr.set_color(text_color)
        new_ctr.move_to(nucleotide_ctr.get_center())

        self.play(
            Transform(nucleotide_ctr, new_ctr),
            TransformFromCopy(first_byte_content[4:6], output_kmer_i[2]),
            run_time=1
        )

        new_ctr = Text("3", font=font, font_size=40, t2w={"weight": BOLD})
        new_ctr.set_color(text_color)
        new_ctr.move_to(nucleotide_ctr.get_center())

        self.play(
            Transform(nucleotide_ctr, new_ctr),
            TransformFromCopy(first_byte_content[6:8], output_kmer_i[3]),
            run_time=1
        )

        self.wait()

        self.play(
            FadeOutToPoint(output_kmer_i, output_kmer_i.get_center() + DOWN),
        )

        output_kmer_i = Text("CGTC", font=font, font_size=60, t2w={"weight": BOLD})
        output_kmer_i.set_color(text_color)
        output_kmer_i.set_stroke(text_color, 2)
        output_kmer_i.next_to(output_kmer, RIGHT, buff=0.5)

        new_ctr = Text("1", font=font, font_size=40, t2w={"weight": BOLD})
        new_ctr.set_color(text_color)
        new_ctr.move_to(nucleotide_ctr.get_center())

        self.wait()

        # Reset the counter
        self.play(
            Transform(nucleotide_ctr, new_ctr)
        )

        self.play(
            TransformFromCopy(first_byte_content[2:4], output_kmer_i[0]),
            run_time = 1
        )

        new_ctr = Text("2", font=font, font_size=40, t2w={"weight": BOLD})
        new_ctr.set_color(text_color)
        new_ctr.move_to(nucleotide_ctr.get_center())

        self.play(
            Transform(nucleotide_ctr, new_ctr),
            TransformFromCopy(first_byte_content[4:6], output_kmer_i[1]),
            run_time=1
        )

        new_ctr = Text("3", font=font, font_size=40, t2w={"weight": BOLD})
        new_ctr.set_color(text_color)
        new_ctr.move_to(nucleotide_ctr.get_center())

        self.play(
            Transform(nucleotide_ctr, new_ctr),
            TransformFromCopy(first_byte_content[6:8], output_kmer_i[2]),
            run_time=1
        )

        new_ctr = Text("0", font=font, font_size=40, t2w={"weight": BOLD})
        new_ctr.set_color(text_color)
        new_ctr.move_to(nucleotide_ctr.get_center())

        self.wait()

        self.play(
            Transform(nucleotide_ctr, new_ctr),
            byte_ptr_group.animate.move_to(second_byte_content[0].get_center() + 3 * DOWN/ 2),
            run_time=1
        )
        self.play(
            TransformFromCopy(second_byte_content[0:2], output_kmer_i[3]),
            run_time=1
        )

        self.wait()

        self.play(
            FadeOutToPoint(output_kmer_i, output_kmer_i.get_center() + DOWN),
        )

        new_ctr = Text("2", font=font, font_size=40, t2w={"weight": BOLD})
        new_ctr.set_color(text_color)
        new_ctr.move_to(nucleotide_ctr.get_center())

        output_kmer_i = Text("GTCA", font=font, font_size=60, t2w={"weight": BOLD})
        output_kmer_i.set_color(text_color)
        output_kmer_i.set_stroke(text_color, 2)
        output_kmer_i.next_to(output_kmer, RIGHT, buff=0.5)

        self.play(
            Transform(nucleotide_ctr, new_ctr),
            byte_ptr_group.animate.move_to(first_byte_content[0].get_center() + 3 * DOWN/ 2),
            run_time=1
        )

        self.play(
            TransformFromCopy(first_byte_content[4:6], output_kmer_i[0]),
            run_time = 1
        )

        new_ctr = Text("3", font=font, font_size=40, t2w={"weight": BOLD})
        new_ctr.set_color(text_color)
        new_ctr.move_to(nucleotide_ctr.get_center())

        self.play(
            Transform(nucleotide_ctr, new_ctr),
            TransformFromCopy(first_byte_content[6:8], output_kmer_i[1]),
            run_time=1
        )


        new_ctr = Text("0", font=font, font_size=40, t2w={"weight": BOLD})
        new_ctr.set_color(text_color)
        new_ctr.move_to(nucleotide_ctr.get_center())
        self.play(
            Transform(nucleotide_ctr, new_ctr),
            byte_ptr_group.animate.move_to(second_byte_content[0].get_center() + 3 * DOWN/ 2),
            run_time=1
        )

        self.play(
            TransformFromCopy(second_byte_content[0:2], output_kmer_i[2]),
            run_time=1
        )

        new_ctr = Text("1", font=font, font_size=40, t2w={"weight": BOLD})
        new_ctr.set_color(text_color)
        new_ctr.move_to(nucleotide_ctr.get_center())

        self.play(
            Transform(nucleotide_ctr, new_ctr),
            TransformFromCopy(second_byte_content[2:4], output_kmer_i[3]),
            run_time=1
        )

        generated_kmers_example = Text("Generated Kmers", font=font, font_size=80, t2w={"weight": BOLD})
        generated_kmers_example.set_color(text_color)
        generated_kmers_example.set_stroke(text_color, 2)

        output_kmers = Text("""
        ACGT
        CGTC
        GTCA
        """, font=font, font_size=60, t2w={"weight": BOLD})
        output_kmers.set_color(text_color)

        generated_kmers_group = VGroup(generated_kmers_example, output_kmers)
        generated_kmers_group.arrange(DOWN, buff=0.5)
        generated_kmers_group.move_to(middle_of_slide)

        self.wait()

        self.play(
            FadeOut(translated_example_group),
            FadeOut(byte_ptr_group),
            FadeOut(info_group),
            FadeOut(output_kmer),
            FadeOut(output_kmer_i),
            TransformFromCopy(dna_example_sequence, generated_kmers_group),
            dna_example_sequence.animate.move_to(middle_of_slide + 2 * UP),
            run_time=1
        )

        self.wait()

        starts_with = Text("Starts with", font=font, font_size=title_text_size, t2w={"weight": BOLD})
        starts_with.set_color(title_text_color)
        starts_with.to_edge(UP)

        self.play(
            TransformMatchingShapes(generate_kmers_header, starts_with),
            FadeOut(generated_kmers_group),
            FadeOut(dna_example_sequence),
            run_time=1
        )

        starts_with_kmer_example = Text(
            """
            Kmer: GATAGA
            """, font=font, font_size=60, t2w={"weight": BOLD}
        )
        starts_with_kmer_example.set_color(text_color)

        starts_with_prefix_example = Text(
            """
            Prefix: GAT
            """, font=font, font_size=60, t2w={"weight": BOLD}
        )
        starts_with_prefix_example.set_color(text_color)

        starts_with_example_group = VGroup(starts_with_kmer_example, starts_with_prefix_example)
        starts_with_example_group.arrange(DOWN, buff=0.5)



        starts_with_kmer_example_bit_representation = Text(
            """
            Kmer: 100011001000
            """, font=font, font_size=60, t2w={"weight": BOLD}
        )
        starts_with_kmer_example_bit_representation.set_color(text_color)

        starts_with_prefix_example_bit_representation = Text(
            """
            Prefix: 100011
            """, font=font, font_size=60, t2w={"weight": BOLD}
        )
        starts_with_prefix_example_bit_representation.set_color(text_color)

        starts_with_example_bit_representation_group = VGroup(starts_with_kmer_example_bit_representation, starts_with_prefix_example_bit_representation)
        starts_with_example_bit_representation_group.arrange(DOWN, buff=0.5)
        starts_with_example_bit_representation_group.move_to(middle_of_slide)
        starts_with_prefix_example_bit_representation.align_to(starts_with_kmer_example_bit_representation, LEFT)

        starts_with_example_group.align_to(starts_with_example_bit_representation_group, LEFT)
        starts_with_example_group.align_to(starts_with_example_bit_representation_group, UP)
        starts_with_prefix_example.align_to(starts_with_kmer_example, LEFT)

        starts_with_bullet_points = Text(
            """
            • Shift the Kmer
            • Compare Kmer with prefix
            """, font=font, font_size=bullet_point_text_size, t2w={"weight": BOLD}
        )
        starts_with_bullet_points.set_color(text_color)
        starts_with_bullet_points.move_to(bullet_point_pos)

        self.play(Write(starts_with_bullet_points))

        self.wait()
        self.play(Write(starts_with_example_group))
        self.wait()

        self.play(
            TransformMatchingShapes(starts_with_example_group, starts_with_example_bit_representation_group),
            run_time=1
        )
        self.wait()
        self.play(
            FadeOut(starts_with_kmer_example_bit_representation[-6:]),
            run_time=1
        )
        self.play(
            starts_with_kmer_example_bit_representation[-12:-6].animate.set_color(GREEN_E),
            starts_with_prefix_example_bit_representation[-6:].animate.set_color(GREEN_E),
            run_time=1
        )

        self.wait()
