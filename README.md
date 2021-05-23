## Gequoteerde functionaliteit

V: Werkend  
-: Deels werkend met gekende problemen (onderaan beschreven)  
X: Niet werkend of niet geïmplementeerd  
(blanco): TODO  


|   | Functionaliteit      | Status |
|---|---------------------------|---|
| 1 | 2D L-systemen             | V |
|   | Met haakjes               | V |
|   | Stochastisch              | X |
| 2 | Transformaties            | V |
|   | Eye-point                 | V |
|   | Projectie                 | V |
| 3 | Platonische Lichamen      | V |
|   | Kegel en cylinder         | V |
|   | Bol                       | V |
|   | Torus                     | V |
|   | 3D L-systemen             | V |
| 4 | Z-buffering (lijnen)      | - |
| 5 | Triangulatie              | V |
|   | Z-buffering (driehoeken)  | V |
| 6 | 3D fractalen              | x |
|   | BuckyBall                 | x |
|   | Mengerspons               | x |
|   | View Frustum              | x |
| 7 | Ambient licht             | V |
|   | Diffuus licht (oneindig)  | V |
|   | Diffuus licht (puntbron)  | - |
|   | Speculair licht           | x |
| 8 | Schaduw                   | x |
|   | Texture mapping           | x |
| 9 | Bollen en cylinders       | x |
|   | UV-coordinaten            | x |
|   | Cube mapping              | x |

Geïmplementeerde vorm van texture mapping: ...

## Gekende problemen 

* ZBuffering met lijnen: in de tekeningen waar enkel kubussen aanwezig zijn, kan men sommige foutjes terug vinden in de overlappingen.
  Dit is mogelijk echter z-fighting in plaats van een fout in de z-buffering met lijnen:
    + z_buffered_wireframes028.ini 
    + z_buffered_wireframes115.ini
  
* Diffuus licht (puntbron) heeft een licht gradient verschil ten opzichte van de gegeven .pngs.
    + diffuus134.ini
    + diffuus217.ini

## Niet-gequoteerde functionaliteit
...

## Extra functionaliteit, niet in de opgaves beschreven
...

