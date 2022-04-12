      COMMON / VTIMTUP / VTXTIM, VTXTER, VTXFTD
      REAL               VTXTIM, VTXTER, VTXFTD

      COMMON / VTIMCOM/ WhatVertex
      INTEGER           WhatVertex

      CHARACTER*(100)  VTIMTUP_form(2)
      DATA             VTIMTUP_form /
     &                 'VTIMTUP',
     &                 'VTXTIM:R,VTXTER:R,VTXFTD:R'/
