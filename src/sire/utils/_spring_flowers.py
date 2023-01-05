from dataclasses import dataclass as _dataclass
from dataclasses import field as _field
from typing import List as _List
from typing import Dict as _Dict

from ._console import Console

__all__ = ["SpringFlowers"]


@_dataclass
class SpringFlowers:
    """This is the colourful 'SpringFlowers' theme"""

    frames: _Dict[int, _List[str]] = _field(default_factory=dict)

    panel_colors: _List[str] = _field(default_factory=list)

    panel_color_count = 0

    def should_highlight(self):
        return False

    def highlighter(self):
        return None

    def should_markup(self):
        return False

    def text(self, style):
        if style == "warning":
            return "magenta"
        elif style == "error":
            return "red"
        elif style == "info":
            return "cyan"
        else:
            return "white"

    def error(self):
        return "red"

    def warning(self):
        return "magenta"

    def info(self):
        return "cyan"

    def spinner_success(self, spinner):
        if Console.supports_emojis():
            spinner.green.ok("✔")
        else:
            spinner.green.ok("Success")

    def spinner_failure(self, spinner):
        if Console.supports_emojis():
            spinner.red.fail("✘")
        else:
            spinner.red.fail("Failed")

    def rule(self, style):
        if style is None:
            return "green"
        elif style == "finish":
            return "magenta"
        elif style == "error":
            return self.error()
        elif style == "warning":
            return self.warning()
        elif style == "info":
            return self.info()
        elif style == "iteration":
            return "cyan"
        else:
            return "cyan"

    def panel_box(self, style):
        from rich import box as _box

        if style == "header":
            return _box.HEAVY_EDGE
        elif style == "command":
            return _box.MINIMAL_HEAVY_HEAD
        else:
            return _box.SQUARE

    def padding_style(self, style):
        if style == "header":
            return "on black"
        elif style == "command":
            return "bold white on #222222"
        else:
            return self.panel(style, advance=False)

    def panel(self, style, advance=True):
        if style is None:
            return "on black"

        elif style == "command":
            return "white on #222222"

        elif style == "alternate":
            if len(self.panel_colors) == 0:
                self.panel_colors = ["blue", "cyan"]
                self.panel_color_count = 0

            color = self.panel_colors[self.panel_color_count]

            if advance:
                self.panel_color_count += 1
                if self.panel_color_count >= len(self.panel_colors):
                    self.panel_color_count = 0

            return f"on {color}"

        elif style == "header":
            return "on purple"

        else:
            return "on black"

    def get_frames(self, width: int = 80):
        """Return the frames used to animate a spinner in a console
        of specified width

        This returns the list of frames plus the timeout between
        the list
        """
        if width in self.frames:
            return self.frames[width]

        frames = []
        frames.append("")

        if Console.supports_emojis():
            bar = "👉 👉 👉 😷  😷  😷 👌 👍 👏 👏 👏 👏 👏 "
        else:
            bar = "-> -> -> #WearAMask :-) :-) :-)  "

        for i in range(1, len(bar), 1):
            frames.append(bar[0:i])

        self.frames[width] = (frames, 50)

        return self.frames[width]
