import h5py, numpy as np
from scipy import sparse, stats
from scipy.sparse.linalg import eigsh
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

layers = ['L1', 'L23', 'L4', 'L5', 'L6']
results = {}

for layer in layers:
    path = f'data/connectivity/nbS1-HEX0-{layer}.h5'
    f = h5py.File(path, 'r')
    indices = f['connectivity/full_matrix/edge_indices/block0_values'][:]
    weights = f['connectivity/full_matrix/edges/block0_values'][:].flatten()
    n       = f['connectivity/full_matrix/vertex_properties/table'].shape[0]
    f.close()

    if n < 10:
        print(f"{layer}: n={n} — muito pequeno, pulando")
        continue

    pre  = indices[:,0]; post = indices[:,1]
    J    = sparse.csr_matrix((weights,(pre,post)), shape=(n,n))
    Js   = (J + J.T) / 2
    s2   = np.var(weights); lp = s2 * 4
    k    = min(30, n-2)
    top  = np.sort(eigsh(Js, k=k, which='LM',
                         return_eigenvectors=False))[::-1]
    out  = top[top > lp]

    # Lei de potência
    if len(out) >= 3:
        rnks = np.arange(1, len(out)+1)
        sl, ic, rv, _, _ = stats.linregress(np.log(rnks), np.log(out))
    else:
        sl, rv = np.nan, np.nan

    # Razão de nível
    gaps = np.diff(top[top>0][::-1]); gaps = gaps[gaps>0]
    r = np.minimum(gaps[:-1]/gaps[1:],
                   gaps[1:]/gaps[:-1]).mean() if len(gaps)>2 else np.nan

    results[layer] = dict(n=n, syn=len(weights),
                          spars=len(weights)/n**2*100,
                          lp=lp, n_out=len(out),
                          lmax=top[0], ratio=top[0]/lp,
                          alpha=sl, r2=rv**2 if not np.isnan(rv) else np.nan,
                          r_goe=r, outliers=out)
    print(f"{layer:4s}  n={n:6d}  syn={len(weights):8d}  "
          f"out={len(out):3d}  α={sl:.3f}  R²={rv**2:.3f}  "
          f"⟨r⟩={r:.3f}  λmax/λ+={top[0]/lp:.1f}×")

# ── Figura: evolução por camada ────────────────────────────────
valid = {k:v for k,v in results.items() if v['n_out'] >= 3}
layer_names = list(valid.keys())

fig, axes = plt.subplots(2, 3, figsize=(15, 9))

# Espectros sobrepostos
ax = axes[0,0]
colors = {'L1':'#888','L23':'#E24B4A','L4':'#EF9F27',
          'L5':'#1D9E75','L6':'#378ADD'}
for ly, res in valid.items():
    out = res['outliers']
    ranks = np.arange(1, len(out)+1)
    ax.plot(ranks, out/out[0], 'o-', color=colors.get(ly,'gray'),
            ms=5, lw=1.5, label=ly)
ax.set_xlabel('Rank'); ax.set_ylabel('λ / λmax (normalizado)')
ax.set_title('Espectros normalizados por camada')
ax.legend()

# Expoente α por camada
ax = axes[0,1]
ly_plot = [l for l in layer_names if not np.isnan(valid[l]['alpha'])]
alphas  = [valid[l]['alpha'] for l in ly_plot]
bars = ax.bar(ly_plot, alphas,
              color=[colors.get(l,'gray') for l in ly_plot],
              alpha=0.8, edgecolor='none')
ax.axhline(-0.416, color='#1D9E75', lw=1.5, ls='--',
           label='MICrONS EM (-0.416)')
ax.axhline(-0.450, color='#E24B4A', lw=1.5, ls=':',
           label='BBP L23 (-0.450)')
ax.set_ylabel('Expoente α')
ax.set_title('Lei de potência por camada')
ax.legend(fontsize=8)
for bar, val in zip(bars, alphas):
    ax.text(bar.get_x()+bar.get_width()/2, val-0.01,
            f'{val:.3f}', ha='center', va='top', fontsize=9)

# ⟨r⟩ por camada
ax = axes[0,2]
rgoe = [valid[l]['r_goe'] for l in ly_plot]
ax.bar(ly_plot, rgoe,
       color=[colors.get(l,'gray') for l in ly_plot],
       alpha=0.8, edgecolor='none')
ax.axhline(0.536, color='gray', lw=1.5, ls='--', label='GOE=0.536')
ax.axhline(0.386, color='gray', lw=1.0, ls=':', label='Poisson=0.386')
ax.axhline(0.654, color='#1D9E75', lw=1.5, ls='--',
           label='MICrONS=0.654')
ax.set_ylabel('⟨r⟩')
ax.set_title('Repulsão de níveis por camada')
ax.legend(fontsize=8)

# N atratores vs N neurônios
ax = axes[1,0]
ns   = [valid[l]['n']     for l in layer_names]
nout = [valid[l]['n_out'] for l in layer_names]
ax.scatter(ns, nout, s=80,
           color=[colors.get(l,'gray') for l in layer_names],
           zorder=5)
for l, x, y in zip(layer_names, ns, nout):
    ax.annotate(l, (x,y), textcoords='offset points',
                xytext=(5,5), fontsize=9)
ax.set_xlabel('N neurônios'); ax.set_ylabel('N atratores')
ax.set_title('N atratores vs tamanho da camada')
sl_n, ic_n, rv_n, _, _ = stats.linregress(np.log(ns), np.log(nout))
x_f = np.linspace(min(ns), max(ns), 100)
ax.plot(x_f, np.exp(ic_n)*x_f**sl_n, '--', color='gray',
        label=f'escala ~ N^{sl_n:.2f}')
ax.legend(fontsize=8)

# λmax/λ+ por camada
ax = axes[1,1]
ratios = [valid[l]['ratio'] for l in layer_names]
ax.bar(layer_names, ratios,
       color=[colors.get(l,'gray') for l in layer_names],
       alpha=0.8, edgecolor='none')
ax.set_ylabel('λmax / λ+')
ax.set_title('Força do atrator dominante')

# Tabela resumo
ax = axes[1,2]; ax.axis('off')
rows = [['Camada','n','Sinapses','N_out','α','R²','⟨r⟩','λmax/λ+']]
for l in layers:
    if l not in results: continue
    r = results[l]
    rows.append([l, f"{r['n']:,}", f"{r['syn']:,}",
                 str(r['n_out']),
                 f"{r['alpha']:.3f}" if not np.isnan(r['alpha']) else '—',
                 f"{r['r2']:.3f}"   if not np.isnan(r['r2'])    else '—',
                 f"{r['r_goe']:.3f}"if not np.isnan(r['r_goe']) else '—',
                 f"{r['ratio']:.1f}×"])
t = ax.table(cellText=rows[1:], colLabels=rows[0],
             loc='center', cellLoc='center')
t.auto_set_font_size(False); t.set_fontsize(8)
t.scale(1, 1.5)
for j in range(len(rows[0])):
    t[0,j].set_facecolor('#2c2c2a')
    t[0,j].set_text_props(color='white', fontweight='bold')
for i, l in enumerate(layers):
    if l in colors:
        for j in range(len(rows[0])):
            t[i+1,j].set_facecolor(colors[l]+'22')

plt.suptitle('TCSV — Análise espectral por camada cortical\n'
             'nbS1-HEX0 (SSCx rato) — BBP Open Data',
             fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('results/layers_spectral.png', dpi=150, bbox_inches='tight')
print("\nSalvo: results/layers_spectral.png")
