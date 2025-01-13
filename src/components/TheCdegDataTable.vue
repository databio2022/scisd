<script setup lang="ts">
import DataTable from 'datatables.net-vue3'
import DataTablesCore from 'datatables.net'
import { ref, defineProps } from 'vue'

DataTable.use(DataTablesCore)

const columns = [
  { data: 'celltype' },
  { data: 'dataset' },
  { data: 'gene' },
  { data: 'avg_log2FC' },
  { data: 'p_val' },
  { data: 'p_val_adj' },
  { data: 'pct_1' },
  { data: 'pct_2' },
  { data: 'final_vs' },
  { data: 'disease' }
]

// 定义 props
const props = defineProps<{
  data: any[]
}>()

const options = ref({
  pageLength: 15,
  lengthMenu: [
    [15, 20, 30, 40, 50, -1],
    ['15', '20', '30', '40', '50', 'All']
  ],
  scrollX: true
})

const layout = {
  title: ''
}
const config = {
  responsive: true
}
</script>

<template>
  <div>
    <DataTable
      :columns="columns"
      :data="props.data"
      class="display"
      :options="options"
      width="100%"
    >
      <thead>
        <tr>
          <th>celltype</th>
          <th>dataset</th>
          <th>gene</th>
          <th>avg_log2FC</th>
          <th>p_val</th>
          <th>p_val_adj</th>
          <th>pct_1</th>
          <th>pct_2</th>
          <th>vs</th>
          <th>disease</th>
        </tr>
      </thead>
      <tfoot>
        <tr>
          <th>celltype</th>
          <th>dataset</th>
          <th>gene</th>
          <th>avg_log2FC</th>
          <th>p_val</th>
          <th>p_val_adj</th>
          <th>pct_1</th>
          <th>pct_2</th>
          <th>vs</th>
          <th>disease</th>
        </tr>
      </tfoot>
    </DataTable>
  </div>
</template>

<style>
@import 'datatables.net-dt';
</style>
