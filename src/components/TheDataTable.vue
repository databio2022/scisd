<script setup lang="ts">
import DataTable from 'datatables.net-vue3'
import DataTablesCore from 'datatables.net'
import { ref, defineProps, computed, onMounted, watch } from 'vue'

import axios from 'axios'

DataTable.use(DataTablesCore)

const columns = [
  { data: 'subject' },
  { data: 'celllabel' },
  { data: 'cluster' },
  { data: 'gene' },
  { data: 'avg_log2FC' },

  { data: 'pct_1' },
  { data: 'pct_2' },

  { data: 'p_val' },
  { data: 'p_val_adj' }
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
  ]
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
          <th>subject</th>
          <th>celllabel</th>
          <th>cluster</th>
          <th>gene</th>
          <th>avg_log2FC</th>
          <th>pct_1</th>
          <th>pct_2</th>

          <th>p_val</th>
          <th>p_val_adj</th>
        </tr>
      </thead>
      <tfoot>
        <tr>
          <th>subject</th>
          <th>celllabel</th>
          <th>cluster</th>
          <th>gene</th>
          <th>avg_log2FC</th>
          <th>pct_1</th>
          <th>pct_2</th>

          <th>p_val</th>
          <th>p_val_adj</th>
        </tr>
      </tfoot>
    </DataTable>
  </div>
</template>

<style>
@import 'datatables.net-dt';
</style>
