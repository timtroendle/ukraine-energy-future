from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw


DARK_GREY = "#424242"
LIGHT_GREY = "#DCDCDC"


def annotate_image(img: Image, zero_line: int) -> Image:
    img = add_margin(img, 100, 0, 0, 0, "white")

    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("Lato-Semibold.ttf", size=32)
    except OSError:
        font = ImageFont.truetype("/cluster/home/trtim/.fonts/Lato-Semibold.ttf", size=32)

    draw.line(
        xy=[(zero_line, 0), (zero_line, 70)],
        fill=LIGHT_GREY,
        width=3
    )

    draw.text(
        xy=(zero_line - 46, 0),
        text="←",
        font=font,
        fill=DARK_GREY,
        align="right"
    )

    draw.text(
        xy=(zero_line + 14, 0),
        text="→",
        font=font,
        fill=DARK_GREY,
        align="left"
    )

    draw.text(
        xy=(zero_line - 410, 0),
        text="Nuclear-and-renewables\n gets cheaper",
        font=font,
        fill=DARK_GREY,
        align="right"
    )

    draw.text(
        xy=(zero_line + 59, 0),
        text="Only-renewables\n gets cheaper",
        font=font,
        fill=DARK_GREY,
        align="left"
    )

    return img


def add_margin(pil_img, top, right, bottom, left, color):
    width, height = pil_img.size
    new_width = width + right + left
    new_height = height + top + bottom
    result = Image.new(pil_img.mode, (new_width, new_height), color)
    result.paste(pil_img, (left, top))
    return result


if __name__ == "__main__":
    img = Image.open(snakemake.input.image)
    img = annotate_image(img, zero_line=snakemake.params.zero_line)
    img.save(snakemake.output.image)
