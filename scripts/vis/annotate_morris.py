from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw


DARK_GREY = "#424242"
LIGHT_GREY = "#DCDCDC"

ZERO_LINE = 916


def annotate_image(img: Image) -> Image:
    img = add_margin(img, 100, 0, 0, 0, "white")

    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype("Lato-Semibold.ttf", size=32)

    draw.line(
        xy=[(ZERO_LINE, 0), (ZERO_LINE, 70)],
        fill=LIGHT_GREY,
        width=3
    )

    draw.text(
        xy=(ZERO_LINE - 46, 0),
        text="←",
        font=font,
        fill=DARK_GREY,
        align="right"
    )

    draw.text(
        xy=(ZERO_LINE + 14, 0),
        text="→",
        font=font,
        fill=DARK_GREY,
        align="left"
    )

    draw.text(
        xy=(ZERO_LINE - 410, 0),
        text="Nuclear-and-renewables\n gets cheaper",
        font=font,
        fill=DARK_GREY,
        align="right"
    )

    draw.text(
        xy=(ZERO_LINE + 59, 0),
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
    img = annotate_image(img)
    img.save(snakemake.output.image)
