import json
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.style.use("ggplot")


def main():
    benchmark_files = sys.argv[1:]
    data = []
    for fpath in benchmark_files:
        with open(fpath) as inhandle:
            results = json.load(inhandle)["results"]
        tool = results[0]["command"].split()[0].split("/")[-1]
        for r in results:
            threads = int(r["parameters"]["threads"])
            for t in r["times"]:
                data.append((tool, threads, float(t)))

    df = pd.DataFrame(data, columns=["tool", "threads", "time"])

    y = "time"
    x = "threads"
    hue = "tool"
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    sns.barplot(data=df, x=x, y=y, hue=hue, ax=ax, dodge=True)
    for b in ax.patches:
        b.set_width(0.39)
    sns.stripplot(data=df, ax=ax, x=x, y=y, hue=hue, dodge=True, color="black", alpha=0.8)
    ax.set(yscale="log")
    yticks = [10, 15, 20, 30, 45, 60, 90, 120, 180, 240]
    ylabels = ["10s", "15s", "20s", "30s", "45s", "1m", "1.5m", "2m", "3m", "4m"]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)

    handles, labels = ax.get_legend_handles_labels()
    del handles[0:2]
    del labels[0:2]
    ax.legend(handles, labels, title="tool")

    tbl = df.groupby(["tool", "threads"])["time"].describe()
    tbl.drop(columns=["count", "25%", "75%"], inplace=True)
    tbl.rename(columns={"50%": "median", "std": "SD"}, inplace=True)
    tbl.reset_index(inplace=True)
    tbl.sort_values(by="threads", inplace=True)
    print(tbl.to_markdown(index=False, floatfmt=".1f"))

    fig.savefig("benchmark.png")


if __name__ == "__main__":
    main()
