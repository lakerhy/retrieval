epsi = linspace(1,9,2);
    for i=2:2
        epsii_snow = 0;
        
      ns = sqrt(epsi(i));

      teta = 0;
      teta = (teta * pi) / 180;
      tei = [asin(sin(teta)/ns);teta]

      [sih,siv] = fresnelc(tei,[epsi(i);1])
    
    end
