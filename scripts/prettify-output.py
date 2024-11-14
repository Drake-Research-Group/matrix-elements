import shutil

def main() -> None:
    terminal_width = shutil.get_terminal_size().columns
    line_char = '\u2500'
    
    line = line_char * terminal_width
    
    print(line)


if __name__ == "__main__":
    main()