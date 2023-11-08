import numpy as np
import math
import tomllib
from loguru import logger
import matplotlib.pyplot as plt

#======== CONSTANTES ========
FICHIER_TOML = "data.toml"
DENSITE_EAU = 997 #[kg/m³]
#============================

#Lecture du fichier toml avec les configurations de la grue
try:
    tomlFile = open(FICHIER_TOML, "rb")
    data = tomllib.load(tomlFile)
    logger.debug("La lecture de data.toml a été effectuée avec succès")
except tomllib.TOMLDecodeError as e:
    logger.error(f"Erreur dans data.toml: {e}")
    exit()

def poussee_archimede(masse):
    return masse * 9.81

def enfoncement(longueur_barge):
    """ 
    Calcule la hauteur de l'enfoncement de la plateforme
    Pré:    masses -> list
            longueur_barge -> int
    Post: retourne un int (hauteur de l'enfoncement)
    """
    return somme_masses() / (DENSITE_EAU * longueur_barge**2)

def somme_masses():
    global data
    sum = 0
    for m in range(len(data["Grue"])):
        sum += int(data["Grue"][str(m)]["masse"])
    sum += data["Barge"]["Masse"]
    return sum

def longueurs_bases(aire_immergee, longueur_barge, angles):
    """
    Calcule la longueur du côté gauche et droit qui est immergée
    Pré:    aire_immergee -> int
            longueur_barge -> int
            angles -> list
    Post: retourne une liste: [longueur_gauche, longueur_droite]
    """
    min_angle = min(angles)
    grande_base = aire_immergee/longueur_barge + (longueur_barge*math.tan(min_angle))/2
    petite_base = aire_immergee/longueur_barge - (longueur_barge*math.tan(min_angle))/2
    return [petite_base, grande_base]

def centre_poussee(longueurs_immergees, longueur_barge, theta):
    lc = (longueur_barge * (longueurs_immergees[1]+2*longueurs_immergees[0]))/(3*(longueurs_immergees[1]+longueurs_immergees[0]))
    hc = (longueurs_immergees[1]**2 + longueurs_immergees[1]*longueurs_immergees[0] + longueurs_immergees[0]**2)/(3*(longueurs_immergees[1]+longueurs_immergees[0]))
    pos_x_prime = -(longueur_barge)/2 + lc
    pos_y_prime = -enfoncement(somme_masses(), longueur_barge) + hc
    pos_x = pos_x_prime*math.cos(theta) - pos_y_prime*math.sin(theta)
    return pos_x

def centre_masse_barge(data):
    CMx = data["Barge"]["Masse"]*((data["Barge"]["Longueur"])/2)
    CMy = data["Barge"]["Masse"]*((data["Barge"]["Hauteur"])/2)
    return (CMx, CMy)

def centre_masse_grue(articulations, angles):
    """
    Calcule le centre de masse dans un plan 2d de la grue dans un axe oxy dont le centre est le début de la base
    Pré: Récupère un dictionnaire sous la forme {'0': {'longueur': 0, 'largeur': 0, 'masse': 0}...}
            et une liste d'angles exprimés en radian contenant n-1 éléments du dictionnaire
    Post: retourne un tuple décrivant la position x,y du centre de masse
    """
    try:
        logger.info("Début du calcul du centre de masse de la grue")
        extr_socle = {"x": articulations["0"]["largeur"]/2, "y": articulations["0"]["longueur"]/2}
        # Initialise les coordonnées des extrémités (avec le socle)
        Ex = [extr_socle["x"]]
        Ey = [extr_socle["y"]]
        # Calcule le x,y de l'extrémité de chaque articulations et le met dans Ex/Ey
        for i in range(1,len(articulations)):
            try:
                Ex.append(articulations[str(i)]["longueur"]*math.cos(angles[i-1]))
                Ey.append(articulations[str(i)]["longueur"]*math.sin(angles[i-1])) 
            except IndexError:
                logger.error(f"Le nombre d'angle n'est pas correct, il faut un total de {len(articulations)-1} angle(s).")
                return
        # Dans ces listes sont stockés les centres de masse
        Cx = [extr_socle["x"]]
        Cy = [extr_socle["y"]]
        # Pour toutes les articulations on calcule le centre de masse
        # (Pas 0 (socle) car déjà calculé dans les tableaux)
        for i in range(1,len(Ex)):
            temp_Cx = 0
            temp_Cy = 0
            # Boucle for pour faire la somme et calculer toutes les articulations
            for u in range(i+1):
                # Si on arrive à l'articulation où l'on cherche le centre de masse
                # on divise l'extrémité par 2 pour obtenir centre de l'articulation
                if i == u:
                    temp_Cx += Ex[u]/2
                    temp_Cy += Ey[u]/2
                    Cx.append(temp_Cx)
                    Cy.append(temp_Cy)
                    break
                # Sinon on additionne simplement les extrémités
                temp_Cx += Ex[u]
                temp_Cy += Ey[u]
        masse_totale = 0
        # Masse * postion pour tous les centres de masse
        pos_massX = 0
        pos_massY = 0
        for i in range(len(Cx)):
            masse_totale += articulations[str(i)]["masse"]
            pos_massX += articulations[str(i)]["masse"] * Cx[i]
            pos_massY += articulations[str(i)]["masse"] * Cy[i]
        # Applique la formule du centre de masse et arrondis à 4 chiffres après la virgule
        CM = (round(pos_massX/masse_totale, 4), round(pos_massY/masse_totale, 4))
        logger.success(f"Le centre de masse de la grue trouvé est {CM}[m]")
        return CM
    except KeyError as e:
        logger.error(f"Vous avez fait une erreur dans le fichier toml, Python le renseigne l'index {e}")
        return

def centre_masse_total(centre_masse_grue,centre_masse_barge, data):
    centre_masse_grue_decale = (centre_masse_grue[0]+data["Barge"]["Declage_grue"]-data["Barge"]["Longueur"]/2, centre_masse_grue[1]+data["Barge"]["Hauteur"]-enfoncement(data["Barge"]["Longueur"]))
    centre_masse_barge_decale = (centre_masse_barge[0]-data["Barge"]["Longueur"]/2, centre_masse_barge[1]-enfoncement(data["Barge"]["Longueur"]))
    masse_articulations_totale = 0
    for i in range(len(data["Grue"])):
        masse_articulations_totale += data["Grue"][str(i)]["masse"]
    CMx = (centre_masse_barge_decale[0] * data["Barge"]["Masse"] + centre_masse_grue_decale[0] * masse_articulations_totale)/(masse_articulations_totale + data["Barge"]["Masse"])
    CMy = (centre_masse_barge_decale[1] * data["Barge"]["Masse"] + centre_masse_grue_decale[1] * masse_articulations_totale)/(masse_articulations_totale + data["Barge"]["Masse"])
    return (CMx, CMy)

def angles_immersion_soulevement(hauteur_enfoncement, hauteur_barge, longueur_barge):
    """
    Calcule les angles d'immersion et de soulèvement
    Pré:    longueur_barge -> int
            hauteur_enfoncement -> int
            hauteur_barge -> int
    Post: retourne une liste avec les deux angles
    """
    angle_immersion = math.atan((2*(hauteur_barge-hauteur_enfoncement))/longueur_barge)
    angle_soulevement = math.atan((2*hauteur_enfoncement)/longueur_barge)
    return [angle_immersion, angle_soulevement]

def aire_immergee(longueur_barge, hauteur_enfoncement):
    """
    Calcule simplement l'aire immergée sous l'eau
    Pré:    longueur_barge -> int
            hauteur_enfoncement -> int
    Post: retourne un int (aire)
    """
    return longueur_barge * hauteur_enfoncement


print(centre_masse_total(centre_masse_grue(data["Grue"], [0,0]), centre_masse_barge(data),data))

xpoints = np.array([0, 7])
ypoints = np.array([3, 3])
plt.plot(xpoints, ypoints, '--')
#plt.show()