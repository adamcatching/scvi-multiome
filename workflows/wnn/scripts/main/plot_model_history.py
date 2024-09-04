import pandas
import matplotlib.pyplot as plt


elbo = pandas.read_csv(snakemake.input.model_history) # type: ignore

val_color = '#D55E00'  # Orange
train_color = '#0072B2'  # Blue

# Plotting
plt.figure(figsize=(8, 6))

# Plot Loss
plt.plot(elbo.index,  elbo['elbo_train'], color=train_color, linestyle='-', label='Training Loss')
plt.plot(elbo.index, elbo['elbo_validation'], color=val_color, linestyle='-', label='Validation Loss')
plt.title('Training and Validation Loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.yscale('log')
plt.legend()
plt.grid(True)

plt.tight_layout()


plt.savefig(snakemake.output.plot_history) # type: ignore