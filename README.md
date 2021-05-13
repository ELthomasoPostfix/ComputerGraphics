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
|   | Z-buffering (driehoeken)  | x |
| 6 | 3D fractalen              |   |
|   | BuckyBall                 |   |
|   | Mengerspons               |   |
|   | View Frustum              |   |
| 7 | Ambient licht             |   |
|   | Diffuus licht (oneindig)  |   |
|   | Diffuus licht (puntbron)  |   |
|   | Speculair licht           |   |
| 8 | Schaduw                   |   |
|   | Texture mapping           |   |
| 9 | Bollen en cylinders       |   |
|   | UV-coordinaten            |   |
|   | Cube mapping              |   |

Geïmplementeerde vorm van texture mapping: ...

## Gekende problemen 

* ZBuffering met lijnen: in de tekeningen waar enkel kubussen aanwezig zijn, kan men sommige foutjes terug vinden in de overlappingen. 
    + z_buffered_wireframes038.ini 
    + z_buffered_wireframes115.ini
* ZBuffering met driehoeken: het crasht met een of andere memory fout wanneer de draw_zbuff_triag() for loop wordt uitgevoerd. Als dit uitgecomment is, dan worden de lege (zwarte achtergrond) .bmp bestanden wel gegenereerd

## Niet-gequoteerde functionaliteit
...

## Extra functionaliteit, niet in de opgaves beschreven
...

