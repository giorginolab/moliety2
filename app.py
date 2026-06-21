from moliety2 import build_demo
from moliety2.ui import APP_CSS, build_theme

demo = build_demo()

if __name__ == "__main__":
    demo.launch(theme=build_theme(), css=APP_CSS)
