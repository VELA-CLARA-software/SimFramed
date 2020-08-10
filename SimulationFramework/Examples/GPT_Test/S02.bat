SET PATH=%PATH%;"C:\Program Files\General Particle Tracer\bin\"
@echo off
rem "gpt.exe" -v -o S02.gdf S02.in GPTLICENSE=***REMOVED***

"gdfa.exe" -o S02_emit.gdf S02_out.gdf position avgx avgy stdx stdBx stdy stdBy stdz stdt nemixrms nemiyrms nemizrms numpar nemirrms avgG avgp stdG avgt avgBx avgBy avgBz CSalphax CSalphay CSbetax CSbetay

"gdfa.exe" -o S02_emit_t.gdf S02_out.gdf time avgx avgy stdx stdBx stdy stdBy stdz stdt nemixrms nemiyrms nemizrms numpar nemirrms avgG avgp stdG avgt avgBx avgBy avgBz CSalphax CSalphay CSbetax CSbetay