function candidato_int = qrs_detection(fm)
%
%candidato_int = qrs_detection(fm)
%
% Función que implementa el algoritmo de Pan y Tompkins para la detección
% de los complejos QRS en una  señal de ECG.
% 
% Entrada: 
%      fm: vector con las fiducial mark.
%
% Salidas:
%
%      candidato_int: vector con las posiciones de los complejos QRS.
%
% J. Pan and W. J. Tompkins, "A real-time QRS detection algorithm". In IEEE
% Trans. on Biomed. Eng., vol 32(3), pp. 230-236. 1985.
%
% 06/11/2014. H. Torres. hmtorres@conicet.gov.ar
% Revisado: 25/10/2019 H. Torres


candidato_int = zeros(size(fm));
noise_peak_int = zeros(size(fm));
n=1;
esperar = 0;
rr_avg_2_int= 0;
rr_missed_limit = 0;

while  n<=length(fm),
  if n==720, % inicialización de los umbrales
      aux = fm(1:n);
      thr_signal_int = max(aux);
      sig_level_int = thr_signal_int;
      aux = aux(aux>0);
      thr_noise_int = min(aux);
      noise_level_int  = thr_noise_int;
  end;

  if (n>720) && (esperar==0),
      if fm(n-1) > thr_signal_int,
          
          candidato_int(n-1) = fm(n-1);
          esperar = 72;
          
          aux = find(candidato_int);
          aux_diff = diff(aux);
          rr_avg_1_int = mean(aux_diff(max(1,end-8):end));
          if rr_avg_2_int== 0,
              rr_avg_2_int = rr_avg_1_int;
          elseif length(aux_diff)>8,
              [aux_n,aux_diff_2] = find(aux_diff > (0.92 * rr_avg_2_int) & aux_diff < (1.16 * rr_avg_2_int));
              if ~isempty(aux_diff_2),
                rr_avg_2_int = mean(aux_diff_2(max(1,end-8):end));
              end;
          end;    
          rr_missed_limit = rr_avg_2_int * 1.66;
          if (length(aux_diff)> 8) & (aux_diff(end) > rr_missed_limit),
              fm_aux = fm(aux(end-1)+1:aux(end)-1);
              aux_fm_aux = aux(end-1)+1 + find(fm_aux > thr_noise_int); 
              if ~isempty(aux_fm_aux),
                  candidato_int(n-1) = 0;
                  candidato_int(aux_fm_aux(1)) = fm(aux_fm_aux(1));
                  n = aux_fm_aux(1) + 1;
                  sig_level_int = 0.25 * fm(n-1) + 0.75 * sig_level_int;
              end;
          else
              sig_level_int = 0.125 * fm(n-1) + 0.875 * sig_level_int;
          end;
          
      else
          noise_level_int = 0.125 * fm(n-1) + 0.875 * noise_level_int;
          noise_peak_int(n-1) = 1;
      end;
      thr_signal_int = noise_level_int + 0.25 * (sig_level_int - noise_level_int);
      thr_noise_int = 0.5 * thr_signal_int;
  elseif  esperar>0, 
      esperar = esperar - 1;
  end;

  n= n+1;
end;





