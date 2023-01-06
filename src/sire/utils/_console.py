from typing import Union as _Union
from typing import List as _List
from typing import IO as _IO

from datetime import datetime as _datetime

from contextlib import contextmanager as _contextmanager


__all__ = ["Console", "Table"]


# Global rich.Console()
_console = None

# Global console theme
_theme = None


class _NullProgress:
    """Null progress to use if user disables progress"""

    def __init__(self):
        pass

    def __enter__(self, *args, **kwargs):
        return self

    def __exit__(self, *args, **kwargs):
        return False

    def add_task(self, *args, **kwargs):
        return 0

    def update(self, *args, **kwargs):
        return


class _Progress:
    def __init__(self, show_limit: int = 50):
        from rich.progress import Progress, BarColumn, TimeRemainingColumn

        self._progress = Progress(
            "[progress.description]{task.description}",
            BarColumn(
                style="bar.back",
                complete_style="bar.complete",
                finished_style="bar.complete",
            ),
            "[progress.percentage]{task.percentage:>3.0f}%",
            TimeRemainingColumn(),
            auto_refresh=False,
        )
        self._last_update = _datetime.now()
        self._show_limit = show_limit
        self._completed = {}
        self._total = {}

        self._have_entered = False
        self._enter_args = None

    def __enter__(self, *args, **kwargs):
        self._enter_args = (args, kwargs)
        return self

    def __exit__(self, *args, **kwargs):
        if self._have_entered:
            self._progress.__exit__(*args, **kwargs)

        return False

    def add_task(self, description: str, total: int):
        if total > self._show_limit:
            if not self._have_entered:
                args, kwargs = self._enter_args
                self._progress.__enter__(*args, **kwargs)
                self._have_entered = True

            task_id = self._progress.add_task(
                description=description, total=total
            )
            self._completed[task_id] = 0
            self._total[task_id] = total
            return task_id
        else:
            return None

    def update(
        self,
        task_id: int,
        completed: float = None,
        advance: float = None,
        force_update: bool = False,
    ):
        if task_id is None:
            return
        elif completed is not None:
            self._completed[task_id] = completed
        elif advance is not None:
            self._completed[task_id] += advance
        else:
            return

        now = _datetime.now()

        if (
            force_update
            or (now - self._last_update).total_seconds() > 0.1
            or self._total[task_id] - self._completed[task_id] < 5
        ):
            self._progress.update(
                task_id=task_id, completed=self._completed[task_id]
            )
            self._progress.refresh()
            self._last_update = now


class _NullSpinner:
    """Null spinner to use if yaspin is not available
    disables spinners"""

    def __init__(self):
        pass

    def __enter__(self, *args, **kwargs):
        return self

    def __exit__(self, *args, **kwargs):
        return False

    def success(self):
        pass

    def failure(self):
        pass


class Table:
    """This a rich.table, with an additional "to_string()" function"""

    def __init__(self, title=None, show_footer=False, show_edge=True):
        from rich.table import Table as _Table

        self._table = _Table(
            title=title, show_edge=show_edge, show_footer=show_footer
        )

    def add_column(
        self, header, justify="center", style=None, no_wrap=False, footer=None
    ):
        """Add a column called 'header', with specified justification,
        style and wrapping
        """
        if footer is not None:
            footer = str(footer)

        self._table.add_column(
            header=header,
            style=style,
            justify=justify,
            no_wrap=no_wrap,
            footer=footer,
        )

    def add_row(self, row):
        """Add the passed row of data to the table"""
        row = [str(x) if x is not None else None for x in row]
        self._table.add_row(*row)

    def to_string(self):
        """Return this table rendered to a string"""
        from rich.console import Console as _Console

        console = _Console()
        output = ""

        for s in console.render(self._table):
            output += s.text

        return output


class Console:
    """This is a singleton class that provides access to printing
    and logging functions to the console. This uses 'rich'
    for rich console printing
    """

    @staticmethod
    def supports_emojis():
        """Return whether or not you can print emojis to this console"""
        import sys

        if sys.platform == "win32":
            return False
        else:
            return True

    @staticmethod
    def set_debugging_enabled(enabled, level=None):
        """Switch on or off debugging output"""
        console = Console._get_console()
        console._debugging_enabled = bool(enabled)
        console._debugging_level = level

    @staticmethod
    def set_theme(theme):
        """Set the theme used for the console - this should be
        one of the themes in metawards.themes
        """
        global _theme

        if isinstance(theme, str):
            if theme.lower().strip() == "simple":
                from ._simple import Simple

                _theme = Simple()
            elif theme.lower().strip() == "default":
                from ._spring_flowers import SpringFlowers

                _theme = SpringFlowers()
        else:
            _theme = theme

    @staticmethod
    def _get_theme():
        global _theme

        if _theme is None:
            from ._spring_flowers import SpringFlowers

            _theme = SpringFlowers()

        return _theme

    @staticmethod
    def _get_console():
        global _console

        if _console is None:
            from rich.console import Console as _Console

            theme = Console._get_theme()

            _console = _Console(
                record=True,
                highlight=theme.should_highlight(),
                highlighter=theme.highlighter(),
                markup=theme.should_markup(),
                log_time=True,
                log_path=True,
                emoji=Console.supports_emojis(),
            )

            _console._use_spinner = True
            _console._use_progress = True
            _console._debugging_enabled = False
            _console._debugging_level = None

            # also install pretty traceback support
            from rich.traceback import install as _install_rich

            _install_rich()

        return _console

    @staticmethod
    @_contextmanager
    def redirect_output(outdir: str, auto_bzip: bool = True):
        """Redirect all output and error to the directory 'outdir'"""
        import os as os
        import bz2
        from rich.console import Console as _Console

        outfile = os.path.join(outdir, "output.txt")

        if auto_bzip:
            outfile += ".bz2"

            OUTFILE = bz2.open(outfile, "wt", encoding="utf-8")
        else:
            OUTFILE = open(outfile, "wt", encoding="utf-8")

        console = Console._get_console()

        if console is None:
            raise AssertionError("The global console should never be None")

        new_out = _Console(
            file=OUTFILE,
            record=False,
            log_time=True,
            log_path=True,
            emoji=Console.supports_emojis(),
        )
        new_out._use_spinner = False
        new_out._use_progress = False
        new_out._debugging_enabled = console._debugging_enabled
        new_out._debugging_level = console._debugging_level
        old_out = console

        global _console
        _console = new_out

        try:
            yield new_out
        finally:
            _console = old_out
            OUTFILE.close()

    @staticmethod
    def debugging_enabled(level: int = None):
        """Return whether debug output is enabled (optionally for
        the specified level) - if not, then
        anything sent to 'debug' (for that level) is not printed
        """
        console = Console._get_console()

        if not console._debugging_enabled:
            return False

        elif console._debugging_level is not None:
            if level is None:
                return True
            else:
                return level <= console._debugging_level
        else:
            return level is None

    @staticmethod
    def _retrieve_name(variable):
        # thanks to scohe001 on stackoverflow
        # https://stackoverflow.com/questions/18425225/getting-the-name-of-a-variable-as-a-string
        import inspect

        callers_local_vars = (
            inspect.currentframe().f_back.f_back.f_locals.items()
        )
        return [
            var_name
            for var_name, var_val in callers_local_vars
            if var_val is variable
        ][-1]

    @staticmethod
    def debug(
        text: str,
        variables: _List[any] = None,
        level: int = None,
        markdown: bool = False,
        **kwargs,
    ):
        """Print a debug string to the console. This will only be
        printed if debugging is enabled. You can also print the
        values of variables by passing them as a list to
        'variables'
        """
        if not Console.debugging_enabled(level=level):
            return

        if hasattr(text, "__call__"):
            # the user passed in a lambda for delayed printing
            text = text()

        if not isinstance(text, str):
            text = str(text)

        if level is not None:
            text = f"Level {level}: {text}"

        console = Console._get_console()

        if markdown:
            from rich.markdown import Markdown as _Markdown

            try:
                text = _Markdown(text)
            except Exception:
                text = _Markdown(str(text))

        console.log(text, justify="center", _stack_offset=2, **kwargs)

        if variables is not None:
            from rich.table import Table
            from rich import box
            import inspect

            # get the local variables in the caller's scope
            callers_local_vars = inspect.currentframe().f_back.f_locals.items()

            table = Table(box=box.MINIMAL_DOUBLE_HEAD, style="on magenta")
            table.add_column(
                "Name", justify="right", style="cyan", no_wrap=True
            )
            table.add_column("Value", justify="left", style="green")

            for variable in variables:
                # get the name of the variable in the caller
                try:
                    # go for the last matching variable in the caller's scope
                    name = [
                        var_name
                        for var_name, var_val in callers_local_vars
                        if var_val is variable
                    ][-1]
                except Exception:
                    name = "variable"

                table.add_row(name, str(variable))

            console.print(table)

    @staticmethod
    def print(
        text: str,
        markdown: bool = False,
        style: str = None,
        markup: bool = None,
        *args,
        **kwargs,
    ):
        """Print to the console"""
        if markdown:
            from rich.markdown import Markdown as _Markdown

            try:
                text = _Markdown(text)
            except Exception:
                text = _Markdown(str(text))

        theme = Console._get_theme()
        style = theme.text(style)

        try:
            Console._get_console().print(text, style=style, markup=markup)
        except UnicodeEncodeError:
            # this output can't cope with a complex theme - switch
            # to theme 'simple'
            Console.set_theme("simple")

            if isinstance(text, str):
                # try to print this as latin-1
                try:
                    text = text.encode("latin-1", errors="replace").decode(
                        "latin-1"
                    )
                    Console._get_console().print(text, markup=markup)
                except Exception:
                    # accept that this isn't going to be printable
                    pass

    @staticmethod
    def rule(title: str = None, style=None, **kwargs):
        """Write a rule across the screen with optional title"""
        from rich.rule import Rule as _Rule

        Console.print("")
        theme = Console._get_theme()
        style = theme.rule(style)
        Console.print(_Rule(title, style=style))

    @staticmethod
    def panel(
        text: str,
        markdown: bool = False,
        width=None,
        padding: bool = True,
        style: str = None,
        expand=True,
        *args,
        **kwargs,
    ):
        """Print within a panel to the console"""
        from rich.panel import Panel as _Panel

        if markdown:
            from rich.markdown import Markdown as _Markdown

            text = _Markdown(text)

        theme = Console._get_theme()
        padding_style = theme.padding_style(style)
        style = theme.panel(style)
        box = theme.panel_box(style)

        from rich.padding import Padding as _Padding

        if padding:
            text = _Padding(text, (1, 2), style=padding_style)
        else:
            text = _Padding(text, (0, 1), style=padding_style)

        Console.print(
            _Panel(
                text,
                box=box,
                width=width,
                expand=expand,
                style=style,
                *args,
                **kwargs,
            )
        )

    @staticmethod
    def error(text: str, *args, **kwargs):
        """Print an error to the console"""
        Console.rule("ERROR", style="error")
        Console.print(text, style="error", *args, **kwargs)
        Console.rule(style="error")

    @staticmethod
    def warning(text: str, *args, **kwargs):
        """Print a warning to the console"""
        Console.rule("WARNING", style="warning")
        Console.print(text, style="warning", *args, **kwargs)
        Console.rule(style="warning")

    @staticmethod
    def info(text: str, *args, **kwargs):
        """Print an info section to the console"""
        Console.rule("INFO", style="info")
        Console.print(text, style="info", *args, **kwargs)
        Console.rule(style="info")

    @staticmethod
    def center(text: str, *args, **kwargs):
        from rich.text import Text as _Text

        Console.print(_Text(str, justify="center"), *args, **kwargs)

    @staticmethod
    def command(text: str, *args, **kwargs):
        Console.panel(text, style="command", padding=False)

    @staticmethod
    def print_profiler(profiler, *args, **kwargs):
        Console.print(str(profiler))

    @staticmethod
    def print_exception():
        """Print the current-handled exception"""
        console = Console._get_console()
        console.print_exception()

    @staticmethod
    def set_use_spinner(use_spinner: bool = True):
        console = Console._get_console()
        console._use_spinner = use_spinner

    @staticmethod
    def set_use_progress(use_progress: bool = True):
        console = Console._get_console()
        console._use_progress = use_progress

    @staticmethod
    def progress(visible: bool = True, show_limit: int = 50):
        console = Console._get_console()

        if visible and console._use_progress:
            return _Progress(show_limit=show_limit)
        else:
            return _NullProgress()

    @staticmethod
    def spinner(text: str = ""):
        try:
            from yaspin import yaspin, Spinner

            have_yaspin = True
        except ImportError:
            have_yaspin = False

        if not have_yaspin:
            return _NullSpinner()

        console = Console._get_console()

        if console._use_spinner:
            theme = Console._get_theme()
            frames, delay = theme.get_frames(
                width=console.width - len(text) - 10
            )
            sp = Spinner(frames, delay)

            y = yaspin(sp, text=text, side="right")

            y.success = lambda: theme.spinner_success(y)
            y.failure = lambda: theme.spinner_failure(y)

            return y
        else:
            return _NullSpinner()

    @staticmethod
    def save(file: _Union[str, _IO]):
        """Save the accumulated printing to the console to 'file'.
        This can be a file or a filehandle. The buffer is
        cleared after saving
        """

        # get the console contents
        try:
            text = Console._get_console().export_text(clear=True, styles=False)
        except Exception as e:
            Console.error(f"Cannot get console output: {e.__class__} {e}")
            return

        if isinstance(file, str):
            with open(file, "w", encoding="UTF-8") as FILE:
                FILE.write(text)
        else:
            try:
                if file.encoding.lower() != "utf-8":
                    # make sure that encoding the file to the right character
                    # set will replace invalid characters, rather than throw
                    # unicode encoding errors
                    text = text.encode(file.encoding, errors="replace").decode(
                        file.encoding
                    )

                file.write(text)
            except UnicodeEncodeError:
                # this still didn't work - use a latin-1 round-robin
                # to make everything ascii...
                text = text.encode("latin-1", errors="replace").decode(
                    "latin-1"
                )
                file.write(text)
