
def convert_to_morse(cyfra):
    #chytra konwersja
    s = ''
    for i in range(0, 5):
        if i < cyfra and i >= cyfra - 5:
            s += '.'
        else:
            s += '-'
    return s

liczby = '0123456789'

def menu():
    czyliczba = False
    liczba = input('Podaj cyfre: ')
    for znak in liczby:
        if str(znak) == str(liczba):
            czyliczba = True
            continue
    if czyliczba:
        morse = convert_to_morse(liczba)
        print("podana liczba w alfabecie morse'a: " + str(morse))
    else:
        print("Podany napis nie jest cyfra!")
        return menu

# to jest wywolanie funkcji menu, tego ci zabraklo
menu()